marker_list <- list(
  "CD4 Naive T cells"           = c("CCR7", "SELL", "TCF7", "LEF1"),
  "CD8 Central Memory T cells"  = c("CD8B","CD248", "NT5E", "NELL2"),
  "CD4 Central Memory T cells"  = c("IL7R", "MAF", "LTB", "GPR183","ANXA1"),
  "CD4 Early Activated T cells" = c("CD69", "FOS", "JUN", "HSPA1B", "EGR1"),
  "CD8 TEMRA T cells"           = c("GZMB", "PRF1", "GNLY", "FCGR3A", "TBX21"),
  "Gamma-delta T cells Vg9Vd2"  = c("TRDV2", "TRGV9", "KLRB1", "NCR3", "KLRG1"),
  "Gamma-delta T cells Vd1"     = c("TRDV1", "TRGC2", "IKZF2", "SOX4", "CXCR3"),
  "Tregs"                       = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2"),
  "ISG-high T cells"            = c("ISG15", "MX1", "IFIT1", "IFIT3", "STAT1")
)
# 按顺序展平基因列表
all_genes <- unique(unlist(marker_list))
# 3. 计算每个细胞类型的平均表达量（按 subcellType 分组）
avg_expr <- AverageExpression(aim_cell, 
                              assays = "RNA", 
                              features = all_genes, 
                              group.by = "subcellType",
                              layer = "data")$RNA

desired_order <- names(marker_list)
common_types <- intersect(desired_order, colnames(avg_expr))
avg_expr <- avg_expr[, common_types, drop = FALSE]
# 对基因进行 z-score 标准化
avg_expr <- as.matrix(avg_expr)
avg_expr_scaled <- t(scale(t(avg_expr)))
avg_expr_scaled[is.na(avg_expr_scaled)] <- 0
# 创建行注释（基因所属细胞类型分组）
gene_to_group <- rep(names(marker_list), times = sapply(marker_list, length))
names(gene_to_group) <- unlist(marker_list)
gene_to_group <- gene_to_group[all_genes]  # 只保留 all_genes 中的基因
annotation_row <- data.frame(Subset = gene_to_group)
rownames(annotation_row) <- all_genes
# 定义分组颜色
subset_colors <- c(
  "CD4 Naive T cells"           = "#1f77b4",
  "CD8 Central Memory T cells"  = "#ff7f0e",
  "CD4 Central Memory T cells"  = "#2ca02c",
  "CD4 Early Activated T cells" = "#d62728",
  "CD8 TEMRA T cells"           = "#9467bd",
  "Gamma-delta T cells Vg9Vd2"  = "#8c564b",
  "Gamma-delta T cells Vd1"     = "#e377c2",
  "Tregs"                       = "#7f7f7f",
  "ISG-high T cells"            = "#bcbd22"
)
annotation_colors <- list(Subset = subset_colors)
# 绘制热图
library(pheatmap)
my_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
pheatmap(avg_expr_scaled,
         color = my_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = annotation_row,
         annotation_names_row = FALSE, 
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         fontsize_col = 10,
         border_color = NA,
         main = "T cell subset marker expression",
         legend = FALSE,  
         filename = "3.2.Tcell_Heatmap.pdf",
         width = 6,
         height = 7
)

###成人
marker_list <- list(
  "CD4 Central Memory T cells"   = c("LTB", "IL7R", "TRADD", "VIM", "PLP2", "CD40LG", "COTL1"),
  "CD4 Early Activated T cells"  = c("AREG", "SOCS3", "CD69", "JUNB", "NR4A1", "SLC2A3"),
  "CD8 TEMRA T cells"            = c("GZMH", "NKG7", "GZMB", "PRF1", "FGFBP2"),
  "CD4 Naive T cells"            = c("CCR7", "SELL", "TCF7", "LEF1", "MYC", "PRKCQ-AS1", "TSHZ2", "CD4"),
  "CD8 Central Memory T cells"   = c("CD8B", "CD8A", "CCR7", "NELL2", "CD248", "ACTN1"),
  "Tregs"                        = c("FOXP3", "IL2RA", "RTKN2", "CTLA4", "TIGIT", "IKZF2"),
  "NK-like T cells"              = c("GNLY", "TYROBP", "KLRD1", "KLRF1", "FCGR3A", "CD244"),
  "Gamma-delta T cells Vg9Vd2"   = c("TRDV2", "TRGV9", "KLRC1", "TRGC1", "KLRB1", "NCR3"),
  "MAIT cells"                   = c("TRAV1-2", "SLC4A10", "KLRB1", "CEBPD", "NCR3"),
  "ISG-high T cells"             = c("ISG15", "IFI6", "HAVCR2", "PDCD1", "CXCR3", "CD27"),
  "CD8 CX3CR1 Effector T cells"  = c("CX3CR1", "ZNF683", "XCL2", "MYOM2", "LINC02384", "FCRL6")
)
# 按顺序展平基因列表
all_genes <- unique(unlist(marker_list))
# 计算每个细胞类型的平均表达量（按 subcellType 分组）
avg_expr <- AverageExpression(aim_cell, 
                              assays = "RNA", 
                              features = all_genes, 
                              group.by = "subcellType",
                              layer = "data")$RNA

desired_order <- names(marker_list)
common_types <- intersect(desired_order, colnames(avg_expr))
avg_expr <- avg_expr[, common_types, drop = FALSE]
# 对基因进行 z-score 标准化
avg_expr <- as.matrix(avg_expr)
avg_expr_scaled <- t(scale(t(avg_expr)))
avg_expr_scaled[is.na(avg_expr_scaled)] <- 0
# 创建行注释（基因所属细胞类型分组）
gene_to_group <- rep(names(marker_list), times = sapply(marker_list, length))
names(gene_to_group) <- unlist(marker_list)
gene_to_group <- gene_to_group[all_genes]  # 只保留 all_genes 中的基因
annotation_row <- data.frame(Subset = gene_to_group)
rownames(annotation_row) <- all_genes
# 定义分组颜色
subset_colors <-  c(
  "CD4 Central Memory T cells"   = "#fde725",
  "CD4 Early Activated T cells"  = "#ff7f0e",
  "CD8 TEMRA T cells"            = "#d62728",
  "CD4 Naive T cells"            = "#1f77b4",
  "CD8 Central Memory T cells"   = "#e08214",
  "Tregs"                        = "#9467bd",
  "NK-like T cells"              = "#bcbd22",
  "Gamma-delta T cells Vg9Vd2"   = "#2ca02c",
  "MAIT cells"                   = "#3b528b",
  "ISG-high T cells"             = "#17becf",
  "CD8 CX3CR1 Effector T cells"  = "#e377c2"
)

annotation_colors <- list(Subset = subset_colors)
# 绘制热图
library(pheatmap)
my_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
pheatmap(avg_expr_scaled,
         color = my_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = annotation_row,
         annotation_names_row = FALSE, 
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         fontsize_col = 10,
         border_color = NA,
         main = "T cell subset marker expression",
         legend = FALSE,  
         filename = "3.2.Tcell_Heatmap.pdf",
         width = 6,
         height = 7
)
