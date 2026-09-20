
library(tidyverse)
library(edgeR)
library(org.Hs.eg.db)

setwd("/GSE279448_bulk")

# ---------- 1. 读取文件并合并表达矩阵 ----------
all_files <- list.files(pattern = "*.htseq.results.gz", full.names = FALSE)
cat("当前路径下共有", length(all_files), "个 htseq 文件\n")

# 提取 donor 编号
donor_from_filename <- all_files %>%
  str_split_fixed("_", n = 2) %>%
  .[, 2] %>%
  str_remove("\\.htseq\\.results\\.gz$")

file_donor_map <- data.frame(filename = all_files, donor = donor_from_filename,
                             stringsAsFactors = FALSE)

# 匹配临床数据中的样本
donors <- clinical$donor
matched <- file_donor_map %>% filter(donor %in% donors)
cat("匹配到", nrow(matched), "个患者样本\n")

# 读取文件函数
read_htseq_file <- function(f, donor_label) {
  read_tsv(f, col_names = c("gene", donor_label), col_types = cols())
}

# 逐个读取并合并
cat("开始读取文件...\n")
expr_list <- list()
for (i in 1:nrow(matched)) {
  df <- read_htseq_file(matched$filename[i], matched$donor[i])
  expr_list[[i]] <- df
  if (i %% 5 == 0) cat("已读取", i, "/", nrow(matched), "\n")
}

bulk_counts <- expr_list %>%
    purrr::reduce(dplyr::full_join, by = "gene") %>%
    dplyr::arrange(gene)

bulk_mat <- bulk_counts %>%
  column_to_rownames("gene") %>%
  as.matrix()
cat("合并后表达矩阵维度：", nrow(bulk_mat), "×", ncol(bulk_mat), "\n")

# ---------- 2. 去除统计行并过滤低表达基因 ----------
valid_genes <- !grepl("^__", rownames(bulk_mat))
bulk_mat <- bulk_mat[valid_genes, ]
cat("去除统计行后，保留", nrow(bulk_mat), "个基因\n")

# 保留至少在 5 个样本中 count >= 5 的基因（可根据需要调整）
keep_genes <- rowSums(bulk_mat >= 5) >= 5
bulk_mat_filtered <- bulk_mat[keep_genes, ]
cat("过滤低表达后，保留", nrow(bulk_mat_filtered), "个基因\n")

# ---------- 3. TMM 标准化 ----------
dge <- DGEList(counts = bulk_mat_filtered)
dge <- calcNormFactors(dge, method = "TMM")
tmm_log2 <- cpm(dge, log = TRUE)
cat("TMM 标准化完成，维度：", dim(tmm_log2), "\n")

# ---------- 4. 基因 ID → Symbol 转换 ----------
gene_ids <- rownames(tmm_log2) %>%
  str_split_fixed("\\.", n = 2) %>%
  .[, 1]

id_to_symbol <- mapIds(org.Hs.eg.db,
                       keys = gene_ids,
                       keytype = "ENSEMBL",
                       column = "SYMBOL",
                       multiVals = "first")

gene_map <- data.frame(
  ensembl_id = gene_ids,
  gene_symbol = id_to_symbol,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(gene_symbol)) %>%
  distinct(ensembl_id, .keep_all = TRUE)

cat("成功转换", nrow(gene_map), "个基因\n")

# ---------- 5. 合并重复 Symbol（取均值）----------
rownames(tmm_log2) <- gene_ids  # 使用纯净的 Ensembl ID
common_ids <- intersect(rownames(tmm_log2), gene_map$ensembl_id)
tmm_subset <- tmm_log2[common_ids, ]

tmm_with_symbol <- as.data.frame(tmm_subset)
tmm_with_symbol$ensembl_id <- rownames(tmm_with_symbol)
tmm_with_symbol <- merge(tmm_with_symbol, gene_map, by = "ensembl_id", all.x = TRUE)

# 按 gene_symbol 分组取均值
cols_to_keep <- !names(tmm_with_symbol) %in% "ensembl_id"
tmm_for_agg <- tmm_with_symbol[, cols_to_keep]
tmm_mean <- aggregate(. ~ gene_symbol, data = tmm_for_agg, FUN = mean, na.rm = TRUE)

rownames(tmm_mean) <- tmm_mean$gene_symbol
tmm_unique <- tmm_mean[, !names(tmm_mean) %in% "gene_symbol"] %>%
  as.matrix()
cat("最终表达矩阵维度（唯一 Symbol）：", dim(tmm_unique), "\n")

# ---------- 6. 关联临床信息并确保列顺序一致 ----------
bulk_clin <- clinical %>%
  filter(donor %in% colnames(tmm_unique)) %>%
  arrange(match(donor, colnames(tmm_unique)))

if (!all(colnames(tmm_unique) == bulk_clin$donor)) {
  tmm_unique <- tmm_unique[, bulk_clin$donor]
}
cat("临床信息匹配：", nrow(bulk_clin), "个样本\n")

# ---------- 7. 保存中间结果 ----------
write.csv(tmm_unique, "2.6.bulk_tmm_log2cpm.csv", row.names = TRUE)
write.csv(bulk_clin, "2.6.bulk_clinical.csv", row.names = FALSE)

# ---------- 8. 计算与 SOFA 的相关性 ----------
bulk_genes <- rownames(tmm_unique)
bulk_assoc <- map_dfr(bulk_genes, function(g) {
  sp <- tryCatch(
    cor.test(as.numeric(tmm_unique[g, ]), bulk_clin$SOFA,
             method = "spearman", exact = FALSE),
    error = function(e) NULL
  )
  if (!is.null(sp)) {
    tibble(gene = g, rho_bulk = sp$estimate, p_bulk = sp$p.value)
  } else {
    NULL
  }
}) %>% arrange(p_bulk)

# FDR 校正
bulk_assoc$fdr_bulk <- p.adjust(bulk_assoc$p_bulk, method = "BH")

# ---------- 9. 筛选集合 B ----------
set_B <- bulk_assoc %>%
  filter(p_bulk < 0.05, abs(rho_bulk) > 0.3)

# ---------- 10. 保存最终结果 ----------
write.csv(set_B, "2.6.bulk.SOFA_gene_filter.csv", row.names = FALSE)


intersect_genes <- intersect(set_A$gene, set_B$gene)
consistent_genes <- intersect_genes[(set_A$rho_partial[match(intersect_genes, set_A$gene)] * set_B$rho_bulk[match(intersect_genes, set_B$gene)]) > 0]

