suppressMessages({
  library(data.table)
  library(jsonlite)
})

# ---------------- 0. 主队列方向字典 (成人模块 consistent genes) ----------------
# 值: +1 = 主队列 SOFA 正相关; -1 = 负相关; IL16/GAS6 为通讯配体 
primary_dir <- c(
  # 成人 T CD8 TEMRA
  UPF3A = -1, SOCS1 = +1, CYB561D1 = -1, PLIN2 = +1, ND4L = +1,
  # 成人髓系 cDC2
  MGA = -1, MTMR9 = -1, SOD1 = +1,
  # 通讯配体 (2026-09-12 起转为正式验证基因)
  IL16 = -1, GAS6 = -1
)
primary_rho <- c(   # 主队列 rho_partial
  UPF3A = -0.347, SOCS1 = 0.351, CYB561D1 = -0.301, PLIN2 = 0.316, ND4L = 0.294,
  MGA = -0.288, MTMR9 = -0.301, SOD1 = 0.248,
  IL16 = -0.174, GAS6 = -0.105
)
exploratory_genes <- character(0)   

target_genes <- c(names(primary_dir), "LINC00211")
stopifnot(length(target_genes) == 11)

# ---------------- 1. 读入数据 ----------------
expr0  <- fread("GSE65682_validation/GSE65682_expr_matrix.csv.gz")
pheno  <- fread("GSE65682_validation/GSE65682_pheno.tsv", sep = "\t")
gene_map <- fromJSON("GSE65682_validation/gene_probe_map.json")

# 表达矩阵: 第一列为探针 ID, 其余列为样本
probe_ids <- expr0[[1]]
samples   <- names(expr0)[-1]
expr_mat  <- as.matrix(expr0[, -1])
rownames(expr_mat) <- probe_ids
stopifnot(all(!duplicated(probe_ids)))

# 表型: gsm 为样本 ID
pheno_gsm <- pheno$gsm

# 脓毒症子集: 有 28 天死亡结局
keep  <- !is.na(pheno$mortality_28)
sub_pheno <- pheno[keep, ]
death <- as.integer(sub_pheno$mortality_28)
cat(sprintf("sepsis subset: n=%d, death=%d, alive=%d\n",
            nrow(sub_pheno), sum(death), nrow(sub_pheno) - sum(death)))

# ---------------- 2. 单基因分析 ----------------
avail_probes <- rownames(expr_mat)
# 仅遍历正式验证基因集 (11 个); 需全部存在于探针映射中
genes <- target_genes[target_genes %in% names(gene_map)]
stopifnot(identical(sort(genes), sort(target_genes)))

mw <- function(x_dead, x_alive) {
  
  w <- wilcox.test(x_dead, x_alive, correct = FALSE)
  U <- unname(w$statistic)
  n1 <- length(x_dead); n2 <- length(x_alive)
  rb <- 1 - 2 * U / (n1 * n2)   # rb>0 = 死亡组更低 (与 Python 一致)
  list(U = U, rb = rb, p = w$p.value)
}

rows <- list()
for (g in genes) {
  ps <- gene_map[[g]]$probes
  ps <- ps[ps %in% avail_probes]
  if (length(ps) == 0) {   
    rows[[length(rows) + 1]] <- data.frame(
      gene = g, n_probes = 0, primary_dir = unname(primary_dir[g]),
      verdict = "not assessable (no probe on platform)", stringsAsFactors = FALSE)
    next
  }
 
  gexpr <- colMeans(expr_mat[ps, , drop = FALSE])
  idx   <- match(sub_pheno$gsm, samples)
  gexpr <- gexpr[idx]
  alive <- gexpr[death == 0]
  dead  <- gexpr[death == 1]
  n1 <- length(dead); n2 <- length(alive)

  w  <- mw(dead, alive)
  z      <- (gexpr - mean(gexpr)) / sd(gexpr)
  delta_z <- mean(z[death == 1]) - mean(z[death == 0])
  fit <- glm(death ~ z, family = binomial)
  or_per_sd <- exp(coef(fit)[2])
  p_logit   <- summary(fit)$coefficients[2, 4]

  d0 <- unname(primary_dir[g])
  direction_match <- (sign(w$rb) == -d0)
  consistent <- direction_match & (w$p < 0.05)
  verdict <- ifelse(consistent, "REPLICATED",
             ifelse(direction_match, "direction consistent (n.s.)", "not consistent"))
  rows[[length(rows) + 1]] <- data.frame(
    gene = g, n_probes = length(ps),
    median_survivors = median(alive), median_non_survivors = median(dead),
    delta_median = median(dead) - median(alive),   
    rank_biserial_r = w$rb, p_wilcoxon = w$p,
    delta_z = delta_z, or_per_sd = or_per_sd, p_logistic = p_logit,
    primary_rho_partial = ifelse(g %in% names(primary_rho), unname(primary_rho[g]), NA),
    primary_dir = d0, exploratory = g %in% exploratory_genes, verdict = verdict,
    stringsAsFactors = FALSE
  )
}
res <- rbindlist(rows, fill = TRUE)

expl   <- !is.na(res$exploratory) & res$exploratory
formal <- !is.na(res$primary_dir) & !expl & !is.na(res$p_wilcoxon)
res$fdr_bh <- NA_real_
if (any(formal)) res$fdr_bh[formal] <- p.adjust(res$p_wilcoxon[formal], method = "BH")

vorder <- c("REPLICATED" = 1, "direction consistent (n.s.)" = 2,
            "not consistent" = 3, "not assessable (no probe on platform)" = 4)
res <- res[order(vorder[res$verdict]), ]

cat("\n===== 单基因验证结果 =====\n")
print(res[, .(gene, n_probes, rank_biserial_r, p_wilcoxon, fdr_bh,
              primary_rho_partial, exploratory, verdict)])

# ---------------- 3. 多基因方向性评分 (加权 z; 正式基因 + 通讯配体) ----------------
score_genes <- c("UPF3A", "SOCS1", "CYB561D1", "PLIN2", "MGA", "MTMR9", "SOD1", "IL16", "GAS6")
score_genes <- score_genes[sapply(score_genes, function(g) {
  any(gene_map[[g]]$probes %in% avail_probes)
})]
Z <- sapply(score_genes, function(g) {
  ps <- gene_map[[g]]$probes; ps <- ps[ps %in% avail_probes]
  colMeans(expr_mat[ps, , drop = FALSE])[match(sub_pheno$gsm, samples)]
})
Zscaled <- scale(Z)                     
score   <- Zscaled %*% primary_dir[score_genes] / sqrt(length(score_genes))
score   <- drop(score)
s_alive <- score[death == 0]; s_dead <- score[death == 1]

w2 <- wilcox.test(s_dead, s_alive, correct = FALSE)   # two-sided
rb2 <- 1 - 2 * unname(w2$statistic) / (length(s_dead) * length(s_alive))
cat(sprintf("\n===== 多基因方向性评分 (n=%d genes) =====\n", length(score_genes)))
cat(sprintf("alive %.3f vs dead %.3f | Mann-Whitney p=%.3g (two-sided) | rb=%.3f\n",
            mean(s_alive), mean(s_dead), w2$p.value, rb2))
cat("score genes:", paste(score_genes, collapse = ", "), "\n")


fwrite(data.frame(gsm = sub_pheno$gsm, death = death, score = score),
       "GSE65682_validation/GSE65682_score_data_R.csv")
cat("saved: GSE65682_validation/GSE65682_score_data_R.csv\n")

# ---------------- 4. 保存 ----------------
fwrite(res, "GSE65682_validation/GSE65682_validation_results_R.csv")
cat("\nsaved: GSE65682_validation/GSE65682_validation_results_R.csv\n")
