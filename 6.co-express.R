#儿童
#髓系细胞亚群
ea_cells <- WhichCells(aim_cell, ident = "Inflammatory Monocytes")
expr_mat <- GetAssayData(aim_cell, slot = "data")[c("WNK1", "CCT7", "HDGF", "FBXO7", "TNIP2", "UBE3C", "CBX1", "TRPV2", "VGLL4", "SUN2", "MVD", "IARS2", "U2AF2", "DCXR", "CSNK2A2", "MPHOSPH6",
  "CCL3","CCL3L1","MIF","IL16","RETN","NAMPT","ANXA1","GRN","LGALS9"), ea_cells, drop=FALSE]
#T细胞亚群
ea_cells <- WhichCells(aim_cell, ident = "CD4 Early Activated T cells")
expr_mat <- GetAssayData(aim_cell, slot = "data")[c("UFD1","IMMT","THOP1","IL16"), ea_cells, drop=FALSE]
#成人
#髓系细胞亚群
ea_cells <- WhichCells(aim_cell, ident = "cDC2")
expr_mat <- GetAssayData(aim_cell, slot = "data")[c("SOD1", "MGA", "MTMR9", "MIF", "IL16", "ANXA1", "GAS6", "GRN", "LGALS9"), ea_cells, drop=FALSE]
#T细胞亚群
ea_cells <- WhichCells(aim_cell, ident = "CD8 TEMRA T cells")
expr_mat <- GetAssayData(aim_cell, slot = "data")[c("UPF3A", "SOCS1", "CYB561D1", "ND4L", "PLIN2", "TNF","MIF", "IL16", "GZMA"), ea_cells, drop=FALSE]

# 自动生成所有两两组合
gene_pairs <- combn(rownames(expr_mat), 2, simplify = FALSE)

# 计算共表达相关函数
calc_corr <- function(g1, g2, mat) {
  cells <- colnames(mat)[mat[g1,] > 0 & mat[g2,] > 0]
  if (length(cells) < 10) return(data.frame(Gene1=g1, Gene2=g2, N=length(cells), Spearman_rho=NA, Spearman_p=NA, Pearson_rho=NA, Pearson_p=NA))
  s <- cor.test(mat[g1,cells], mat[g2,cells], method="spearman")
  p <- cor.test(mat[g1,cells], mat[g2,cells], method="pearson")
  data.frame(Gene1=g1, Gene2=g2, N=length(cells), Spearman_rho=s$estimate, Spearman_p=s$p.value, Pearson_rho=p$estimate, Pearson_p=p$p.value)
}

res <- do.call(rbind, lapply(gene_pairs, function(pair) calc_corr(pair[1], pair[2], expr_mat)))

write.csv(res, "3.8.coexpression.T.csv", row.names=FALSE)
print(res)

res_filtered <- res[!is.na(res$Spearman_rho) & !is.na(res$Pearson_rho) & 
                    res$Spearman_rho > 0.5 & res$Pearson_rho > 0.5 & 
                    res$Spearman_p < 0.05 & res$Pearson_p < 0.05, ]
write.csv(res_filtered, "2.8.coexpression_filtered.M.csv", row.names=FALSE)


