
suppressMessages({
  library(data.table)
  library(jsonlite)
})

# ---------------- 0. 主队列方向字典 (儿童全部模块基因) ----------------
# 值: +1 = 主队列 SOFA 正相关 (预期死亡组更高, rb<0); 0 = 探索性; NA = 无探针不可评估
MYELOID <- c("WNK1","CCT7","HDGF","FBXO7","TNIP2","UBE3C","CBX1","TRPV2","VGLL4",
             "SUN2","MVD","IARS2","U2AF2","DCXR","CSNK2A2","MPHOSPH6")
TCELL   <- c("UFD1","IMMT","THOP1")

primary_dir <- c(
  setNames(rep(1L, length(MYELOID)), MYELOID),   # 16 儿童髓系一致基因
  setNames(rep(1L, length(TCELL)),   TCELL),     # 3 儿童T 一致基因
  IL16 = 1L,                                     # 跨谱系 hub: 与全部儿童模块基因正耦合, 继承 dir=+1
  LINC00211 = NA_integer_                        # lncRNA, 平台无探针
)
primary_rho <- c(
  WNK1 = 0.785, CCT7 = 0.679, HDGF = 0.712, FBXO7 = 0.653, TNIP2 = 0.526,
  UBE3C = 0.732, CBX1 = 0.735, TRPV2 = 0.576, VGLL4 = 0.641, SUN2 = 0.562,
  MVD = 0.568, IARS2 = 0.544, U2AF2 = 0.626, DCXR = 0.500, CSNK2A2 = 0.750,
  MPHOSPH6 = 0.579,
  UFD1 = 0.583, IMMT = 0.681, THOP1 = 0.564,
  IL16 = 0.359
)
module <- c(
  setNames(rep("Pediatric myeloid (Inflammatory Monocytes)", length(MYELOID)), MYELOID),
  setNames(rep("Pediatric T (CD4 Early Activated)", length(TCELL)), TCELL),
  IL16 = "Cross-lineage hub",
  LINC00211 = "Myeloid-specific lncRNA"
)

# ---------------- 1. 读入 GSE66099 统一 gcRMA 表达矩阵 ----------------
f_gz <- "GSE66099_validation/GSE66099_series_matrix.txt.gz"
con  <- gzfile(f_gz, "rt", encoding = "UTF-8")
ln   <- 0; begin <- NULL
repeat {
  chunk <- readLines(con, n = 200)
  if (length(chunk) == 0) break
  hit <- grep("!series_matrix_table_begin", chunk, fixed = TRUE)
  if (length(hit)) { begin <- ln + hit[1]; break }
  ln <- ln + length(chunk)
}
close(con)
stopifnot(!is.null(begin))
dt  <- fread(f_gz, skip = begin, sep = "\t", header = TRUE,
             check.names = FALSE)          
probe_ids <- dt[[1]]; samples <- names(dt)[-1]
expr_mat  <- as.matrix(dt[, -1]); rownames(expr_mat) <- probe_ids
cat(sprintf("GSE66099 统一矩阵: %d 探针 x %d 样本\n", nrow(expr_mat), ncol(expr_mat)))

# ---------------- 2. 结局匹配 (Reanalysis of 回溯) ----------------
gsm_map <- fromJSON("GSE66099_validation/gsm_to_original.json")
out_map <- fromJSON("GSE66099_validation/pediatric_outcome_map.json")
gse_gsm <- gsm_map$gsm_order
disease <- unlist(gsm_map$disease)
original<- unlist(gsm_map$original)

shock <- gse_gsm[disease[gse_gsm] == "SepticShock" & gse_gsm %in% names(original)]
cat(sprintf("GSE66099 SepticShock 样本: %d\n", length(shock)))

outcome <- list()
for (g in shock) {
  o <- original[[g]]
  for (gse in names(out_map)) {          
    if (o %in% names(out_map[[gse]])) {
      out <- out_map[[gse]][[o]]
      if (out %in% c("Survivor", "Nonsurvivor")) {
        outcome[[g]] <- ifelse(out == "Nonsurvivor", 1L, 0L)
        break
      }
    }
  }
}
use_gsm <- names(outcome)
death   <- unlist(outcome)
cat(sprintf("匹配到结局的样本: %d (死亡 %d / 存活 %d)\n",
            length(use_gsm), sum(death), length(use_gsm) - sum(death)))

# ---------------- 3. 单基因分析 (多探针均值, GSE66099 自身矩阵) ----------------
gene_map <- fromJSON("GSE66099_validation/gene_probe_map_ped.json")
targets  <- names(gene_map)
avail    <- rownames(expr_mat)
col_idx  <- match(use_gsm, samples)

sub <- data.frame(death = death, row.names = use_gsm)
no_probe <- character(0)
for (g in targets) {
  ps <- gene_map[[g]]$probes
  ps <- ps[ps %in% avail]
  if (length(ps) == 0) {                 
    no_probe <- c(no_probe, g)
    next
  }
  sub[[g]] <- colMeans(expr_mat[ps, col_idx, drop = FALSE])
}
if (length(no_probe)) cat(sprintf("无探针不可评估: %s\n", paste(no_probe, collapse = ", ")))

mw <- function(x_dead, x_alive) {
  w <- wilcox.test(x_dead, x_alive, correct = FALSE)
  U <- unname(w$statistic); n1 <- length(x_dead); n2 <- length(x_alive)
  list(U = U, rb = 1 - 2 * U / (n1 * n2), p = w$p.value)
}

rows <- list()
for (g in targets) {
  d0 <- unname(primary_dir[g])
  if (!(g %in% names(sub))) {            
    rows[[length(rows) + 1]] <- data.frame(
      gene = g, module = unname(module[g]), n_probes = 0L,
      primary_dir = d0, verdict = "not assessable (no probe on platform)",
      stringsAsFactors = FALSE
    )
    next
  }
  gexpr <- sub[[g]]
  alive <- gexpr[sub$death == 0]; dead <- gexpr[sub$death == 1]
  w <- mw(dead, alive)
  z <- (gexpr - mean(gexpr)) / sd(gexpr)
  delta_z <- mean(z[sub$death == 1]) - mean(z[sub$death == 0])
  fit <- glm(sub$death ~ z, family = binomial)
  or_per_sd <- exp(coef(fit)[2])
  n_probes <- sum(gene_map[[g]]$probes %in% avail)
  if (is.na(d0)) {
    verdict <- "not assessable (no probe on platform)"
  } else if (d0 == 0) {
    verdict <- "exploratory"
  } else {
    direction_match <- (sign(w$rb) == -d0)   
    consistent <- direction_match & (w$p < 0.05)
    verdict <- ifelse(consistent, "REPLICATED",
               ifelse(direction_match, "direction consistent (n.s.)", "not consistent"))
  }
  rows[[length(rows) + 1]] <- data.frame(
    gene = g, module = unname(module[g]), n_probes = n_probes,
    median_survivors = median(alive), median_non_survivors = median(dead),
    delta_median = median(dead) - median(alive),   # = non-survivors - survivors
    rank_biserial_r = w$rb, p_wilcoxon = w$p,
    delta_z = delta_z, or_per_sd = or_per_sd,
    primary_rho_partial = unname(primary_rho[g]),
    primary_dir = d0, verdict = verdict,
    stringsAsFactors = FALSE
  )
}
res <- rbindlist(rows, fill = TRUE)

formal <- !is.na(res$primary_dir) & res$primary_dir != 0
res$fdr_bh <- NA_real_
if (any(formal)) res$fdr_bh[formal] <- p.adjust(res$p_wilcoxon[formal], method = "BH")

vorder <- c("REPLICATED" = 1, "direction consistent (n.s.)" = 2,
            "exploratory" = 3, "not consistent" = 4,
            "not assessable (no probe on platform)" = 5)
res <- res[order(vorder[res$verdict]), ]

cat("\n===== 单基因验证结果 =====\n")
print(res[, .(gene, module, n_probes, rank_biserial_r, p_wilcoxon, fdr_bh,
              primary_rho_partial, verdict)])

# ---------------- 4. 多基因方向性评分 (全部一致基因, dir != 0) ----------------
score_genes <- names(primary_dir)[!is.na(primary_dir) & primary_dir != 0]
Zscaled <- scale(sub[, score_genes, drop = FALSE])
score   <- drop(Zscaled %*% primary_dir[score_genes] / sqrt(length(score_genes)))
s_alive <- score[sub$death == 0]; s_dead <- score[sub$death == 1]
w2 <- wilcox.test(s_dead, s_alive, correct = FALSE)   # two-sided
cat(sprintf("\n===== 多基因方向性评分 (n=%d genes, dir=+1) =====\n", length(score_genes)))
cat(sprintf("alive %.3f (SD %.3f) vs dead %.3f (SD %.3f) | Mann-Whitney p=%.3g\n",
            mean(s_alive), sd(s_alive), mean(s_dead), sd(s_dead), w2$p.value))

fwrite(data.frame(gsm = use_gsm, death = death, score = score),
       "GSE66099_validation/GSE66099_score_data_R.csv")
cat("saved: GSE66099_validation/GSE66099_score_data_R.csv\n")

# ---------------- 5. 保存 ----------------
fwrite(res, "GSE66099_validation/GSE66099_validation_results_R.csv")
cat("\nsaved: GSE66099_validation/GSE66099_validation_results_R.csv\n")
