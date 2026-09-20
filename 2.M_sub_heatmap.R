marker_list <- list(
  "Tissue-resident Macrophages"          = c("C1QA", "C1QB", "C1QC", "MS4A4A", "FCGR3A", "CD68"),
  "Neutrophils"                          = c("S100A8", "S100A9", "RETN", "PADI4", "HP", "ALOX5AP"),
  "Inflammatory Monocytes"               = c("EREG", "ID1", "CXCL8", "CCL3", "FOS", "JUN"),
  "IFN-responsive Monocytes"             = c("IFI6", "MX1", "IRF7", "BST2","RSAD2"),
  "Proliferating Neutrophil Progenitors" = c("TOP2A", "MCM5", "PCNA", "STMN1", "MKI67"),
  "Stressed Monocytes"                   = c("TCF7L2","MALAT1","CD163", "TNFRSF1B", "NAMPT", "TXNIP"),
  "cDC2"                                 = c("CD1C", "CLEC10A", "FLT3", "HLA-DRA", "RGS1")
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
subset_colors <- c(
  "Tissue-resident Macrophages"          = "#440154",   
  "Neutrophils"                          = "#cc4778",   
  "Inflammatory Monocytes"               = "#21918c",   
  "IFN-responsive Monocytes"             = "#ff7f0e",   
  "Proliferating Neutrophil Progenitors" = "#3b528b", 
  "Stressed Monocytes"                   = "#5ec962",   
  "cDC2"                                 = "#9467bd"  
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
         main = "Myeloid cell subset marker expression",
         legend = FALSE,  
         filename = "2.2.Mcell_Heatmap.pdf",
         width = 6,
         height = 7
)
###成人
marker_list <- list(
  "Classical Monocytes"                  = c("RETN", "MNDA", "PLBD1","CDA", "CSTA", "MSRB1", "CYP1B1"),
  "Inflammatory Monocytes"               = c("IL1B", "CXCL8", "CCL3", "CCL3L1", "CXCL2", "EREG", "G0S2"),
  "Non-classical Monocytes"              = c("FCGR3A", "MS4A7", "CX3CR1", "CDKN1C", "PELATON", "LST1", "AIF1", "CSF1R"),
  "IFN-responsive Monocytes"             = c("ISG15", "IFI6", "MX1", "STAT1", "LY6E", "OAS1"),
  "Activated Monocytes"                  = c("EGR1", "GADD45B", "GADD45G", "IER2", "DNAJB1", "HSPA1A"),
  "Tissue-resident Macrophages"          = c("C1QA", "C1QB", "C1QC", "MS4A4A", "FCGR3A", "CD68"),
  "cDC2"                                 = c("CD1C", "CLEC10A", "FCGR2B", "GPR183", "FLT3", "HLA-DRA"),
  "Proliferating Neutrophil Progenitors" = c("MKI67", "PTTG1", "STMN1", "CENPF", "BIRC5", "HMGN2"),
  "Neutrophils"                          = c("MMP9", "CD177", "FCGR3B", "CSF3R", "CXCR2", "PADI4"),
  "Stressed Monocytes"                   = c("PET100", "VCAN", "CLU", "MCL1", "COX3", "ND2")
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
subset_colors <- c(
  "Classical Monocytes"                  = "#fde725",
  "Inflammatory Monocytes"               = "#21918c",
  "Non-classical Monocytes"              = "#1f77b4",
  "IFN-responsive Monocytes"             = "#ff7f0e",
  "Activated Monocytes"                  = "#e377c2",
  "Tissue-resident Macrophages"          = "#440154",
  "cDC2"                                 = "#9467bd",
  "Proliferating Neutrophil Progenitors" = "#3b528b",
  "Neutrophils"                          = "#cc4778",
  "Stressed Monocytes"                   = "#5ec962"
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
         main = "Myeloid cell subset marker expression",
         legend = FALSE,  
         filename = "2.2.Mcell_Heatmap.pdf",
         width = 6,
         height = 7
)
