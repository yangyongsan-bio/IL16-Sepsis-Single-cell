aim_cell=readRDS("sub.adult.T.celltype.Rds")

target_clusters <- "CD8 TEMRA T cells"

C_cells <- subset(aim_cell, subset = subcellType == target_clusters)

# 按 donor 计算平均表达（log-normalized data）
avg_exp_list <- AverageExpression(C_cells, assays = "RNA", group.by = "orig.ident",
                                  slot = "data", verbose = FALSE)
avg_exp <- avg_exp_list$RNA  # 提取 RNA 矩阵（基因 × donor）

# 转置为 donor × gene 矩阵
avg_exp_dense <- as.matrix(avg_exp)
avg_exp_t <- t(avg_exp_dense) %>% as.data.frame()
avg_exp_t$donor <- rownames(avg_exp_t)

# 合并临床信息（SOFA）
dat_gene <- left_join(avg_exp_t, clinical %>% dplyr::select(donor, SOFA, Age, inf_group), 
                      by = "donor")

# ---------- 2. 基因预过滤 ----------
# 条件1：至少在 一半的 donor 中表达（>0）
n_samples <- ncol(avg_exp)
min_donors <- floor(n_samples / 2) + 1
n_donors_expr <- rowSums(avg_exp > 0)
# 条件2：表达量有足够的变异（CV > 0.1）
gene_cv <- apply(avg_exp, 1, function(x) sd(x) / mean(x))
# 条件3：平均表达 > 0.1（避免极低表达基因）
gene_mean <- rowMeans(avg_exp)

genes_pass_filter <- names(n_donors_expr[n_donors_expr >= min_donors & 
                                          gene_mean > 0.1 & 
                                          gene_cv > 0.1])

cat("通过预过滤的基因数：", length(genes_pass_filter), "\n")

# ---------- 3. 逐基因关联 SOFA（Spearman，加混杂） ----------
# 主要用 Spearman 单变量，混杂作为敏感性分析
gene_assoc <- map_dfr(genes_pass_filter, function(g) {
  # 单变量 Spearman
  sp <- tryCatch(
    cor.test(dat_gene[[g]], dat_gene$SOFA, method = "spearman", exact = FALSE),
    error = function(e) NULL
  )
  
  # 加混杂的偏相关（Age + inf_group）
  partial_sp <- tryCatch({
    # 使用 ppcor 包的 pcor 或手动残差法
    # 简化：先对 Age 和 inf_group 回归取残差
    res_g <- residuals(lm(as.formula(paste(g, "~ Age + inf_group")), data = dat_gene))
    res_sofa <- residuals(lm(SOFA ~ Age + inf_group, data = dat_gene))
    cor.test(res_g, res_sofa, method = "spearman", exact = FALSE)
  }, error = function(e) NULL)
  
  tibble(
    gene = g,
    rho_raw = ifelse(is.null(sp), NA, sp$estimate),
    p_raw = ifelse(is.null(sp), NA, sp$p.value),
    rho_partial = ifelse(is.null(partial_sp), NA, partial_sp$estimate),
    p_partial = ifelse(is.null(partial_sp), NA, partial_sp$p.value)
  )
})

# FDR 校正
gene_assoc$fdr_raw <- p.adjust(gene_assoc$p_raw, method = "BH")
gene_assoc$fdr_partial <- p.adjust(gene_assoc$p_partial, method = "BH")

# 按原始 p 值排序
gene_assoc <- gene_assoc %>% arrange(p_raw)

write.csv(gene_assoc, file="3.5.T.SOFA_gene.csv", row.names =F,quote=F)

set_A <- gene_assoc %>%
  filter(p_raw < 0.05,                     # 原始相关显著
         abs(rho_raw) > 0.3,               # 效应量足够大
         !is.na(p_partial),                # 偏相关可计算
         p_partial < 0.05,                 # 偏相关显著（控制混杂后仍显著）
         sign(rho_raw) == sign(rho_partial))  # 方向一致
write.csv(set_A, file="3.5.T.SOFA_gene_filter.csv", row.names =F,quote=F)