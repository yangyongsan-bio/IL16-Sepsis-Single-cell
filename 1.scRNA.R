library(Seurat)
library(scDblFinder)
library(dplyr)
library(ggplot2)
library(clustree)
library(harmony)
library(Matrix)
library(org.Hs.eg.db)
library(tidyverse)
library(broom)
setwd("/home/sc")
dir=dir()
scRNAlist <- list()
for (i in seq_along(dir)) {
  cat("正在处理第", i, "个样本:", dir[i], "\n")
  # 1. 读取10X数据（基因名来自features.tsv第一列，即Ensembl ID）
  counts <- Read10X(data.dir = dir[i], gene.column = 1)
  # 2. 创建Seurat对象
  obj <- CreateSeuratObject(counts = counts, project = dir[i])
  # 3. 添加样本标识并重命名细胞
  obj$orig.ident <- dir[i]
  new_cell_names <- paste0(dir[i], colnames(obj))
  obj <- RenameCells(object = obj, new.names = new_cell_names)
  # 4. 基因名转换：Ensembl ID → Gene Symbol
  current_genes <- rownames(obj)
  # 去掉版本号后缀（如 .1, .2）
  clean_ids <- sub("\\.[0-9]+$", "", current_genes)
  # 查询基因符号，未匹配的保留原ID
  gene_symbol <- mapIds(org.Hs.eg.db,
                        keys = clean_ids,
                        keytype = "ENSEMBL",
                        column = "SYMBOL",
                        multiVals = "first")
  gene_symbol[is.na(gene_symbol)] <- current_genes[is.na(gene_symbol)]
  # 5. 替换行名（注意：不能直接修改Seurat对象的行名，需要提取矩阵再重建）
  raw_counts <- GetAssayData(obj, slot = "counts")
  rownames(raw_counts) <- gene_symbol
  # 6. 处理重复基因名（多个Ensembl ID映射到同一Symbol）
  dup <- duplicated(rownames(raw_counts))
  if (any(dup)) {
    cat("  发现重复基因名，正在聚合...\n")
    # 将稀疏矩阵转为稠密矩阵（单个样本矩阵不大，安全）
    dense_mat <- as.matrix(raw_counts)
    # 按基因名分组求和（也可改为求均值，这里用求和）
    raw_counts <- rowsum(dense_mat, group = rownames(dense_mat), reorder = FALSE)
    # 转回稀疏矩阵节省内存
    raw_counts <- as(raw_counts, "dgCMatrix")
  }
  # 7. 用转换后的矩阵重建Seurat对象（保留原有meta.data）
  obj <- CreateSeuratObject(counts = raw_counts,
                            project = dir[i],
                            meta.data = obj@meta.data)
  # 8. 存入列表
  scRNAlist[[i]] <- obj
  cat("  完成，细胞数:", ncol(obj), "基因数:", nrow(obj), "\n")
}
scRNA <- Reduce(function(x, y) merge(x, y), scRNAlist[1:133])
saveRDS(scRNA, file="total.Rds")

setwd("/home/Seurat")
library(readxl)
library(dplyr)
# 1. 读取样本表（两个sheet）
adult <- read_excel("Subsample.xlsx", sheet = 1)
pediatric <- read_excel("Subsample.xlsx", sheet = 2)
# 2. 提取成人Sepsis和儿童Sepsis的Donor ID列表（直接取整列）
adult_ids <- adult$`Donor.ID`
ped_ids <- pediatric$`Donor.ID`
# 3. 从scRNA中subset出成人Sepsis细胞
scRNA$orig.ident[scRNA$orig.ident == "Skin-152"] <- "Skin-S152" #有个样本名错误
scRNA_adult <- subset(scRNA, subset = orig.ident %in% adult_ids)
# 4. 提取成人样本表的临床信息
# 选择需要添加的列：取 adult的所有列名，然后去掉 "Donor ID"这一列
cols_adult <- setdiff(names(adult), "Donor.ID")
#  使用 left_join 将临床信息添加到成人子集的meta.data
scRNA_adult@meta.data <- scRNA_adult@meta.data %>%
  left_join(adult[, c("Donor.ID", cols_adult)], by = c("orig.ident" = "Donor.ID"))
rownames(scRNA_adult@meta.data)=colnames(scRNA_adult)
saveRDS(scRNA_adult,"adult.Rds")
# 同理儿童也一样
scRNA_pediatric <- subset(scRNA, subset = orig.ident %in% ped_ids)
cols_pediatric <- setdiff(names(pediatric), "Donor.ID")
scRNA_pediatric@meta.data <- scRNA_pediatric@meta.data %>%
  left_join(pediatric[, c("Donor.ID", cols_pediatric)], by = c("orig.ident" = "Donor.ID"))
rownames(scRNA_pediatric@meta.data)=colnames(scRNA_pediatric)
saveRDS(scRNA_pediatric,"pediatric.Rds")

############### scRNA-seq分析 ###############
scRNA=readRDS("pediatric.Rds")
scRNA = as.SingleCellExperiment(scRNA)
scRNA <- scDblFinder(scRNA, samples="orig.ident")
scRNA=as.Seurat(scRNA)
####质控####
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
pdf("p1.pdf")
scRNA@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill=orig.ident)) + 
  geom_density(alpha = 0.15) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("log10 Counts density") +
  geom_vline(xintercept = 100) +
  geom_vline(xintercept = 30000)
dev.off()
pdf("p2.pdf")
scRNA@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() +
  theme_classic() +
  ylab("log10 gene density") +
  geom_vline(xintercept = 200)+
  geom_vline(xintercept = 6000)+
  geom_vline(xintercept = 10000)
dev.off()
pdf("p3.pdf")
scRNA@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  geom_vline(xintercept = 5) +
  geom_vline(xintercept = 10)
dev.off()
#儿童的
scRNA = subset(scRNA,cells=which(scRNA$nFeature_RNA > 200 & scRNA$nFeature_RNA < 6000 & scRNA$percent.mt < 5 
& scRNA$nCount_RNA > 500 & scRNA$nCount_RNA < 30000 & scRNA$scDblFinder.class=="singlet"))
#成人的
scRNA = subset(scRNA,cells=which(scRNA$nFeature_RNA > 100 & scRNA$nFeature_RNA < 10000 & scRNA$percent.mt < 5 
& scRNA$nCount_RNA > 500 & scRNA$nCount_RNA < 100000 & scRNA$scDblFinder.class=="singlet"))

saveRDS(scRNA, "adult.qc.Rds")
####标准化数据####
scRNA <- NormalizeData(scRNA, normalization.method = 'LogNormalize', scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA <- ScaleData(scRNA, vars.to.regress = "percent.mt", features = VariableFeatures(scRNA))
scRNA <- RunPCA(scRNA)
####去除批次####
#####直接运行harmony报错。
## 提取PCA嵌入矩阵
pca_embedding <- Embeddings(scRNA, reduction = "pca")
## 运行Harmony而不依赖Seurat对象
harmony_embedding <- HarmonyMatrix(
  data_mat = pca_embedding,
  meta_data = scRNA@meta.data,
  vars_use = "orig.ident",
  do_pca = FALSE  # 因为我们已经提供了PCA嵌入
)
## 将结果添加回Seurat对象
scRNA[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_embedding,
  key = "harmony_",
  assay = "RNA"
)
ggsave('ElbowPlot.pdf', ElbowPlot(scRNA,ndims = 50), width = 10, height = 5)
pc.num=1:30
scRNA <- RunUMAP(scRNA,reduction="harmony", dims=pc.num)
scRNA <- RunTSNE(scRNA,reduction="harmony", dims=pc.num)
scRNA <- FindNeighbors(scRNA,reduction="harmony",dims = pc.num,k.param = 20)
scRNA <- FindClusters(scRNA, resolution = seq(0.1,0.5,by=0.1))
ggsave('clustree.pdf', clustree(scRNA), width = 12, height = 10)
scRNA$Seurat_clusters=scRNA$RNA_snn_res.0.5
Idents(scRNA)=scRNA$Seurat_clusters
saveRDS(scRNA, "adult.0.5.Rds")
allMarkers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(allMarkers, file="adult.allMarkers.csv", row.names =F,quote=F)

###成人注释
scRNA=readRDS("adult.0.5.Rds")
scRNA=subset(scRNA,cells=which(scRNA$seurat_clusters != 13)) #红细胞
scRNA=subset(scRNA,cells=which(scRNA$seurat_clusters != 15)) #混合细胞
scRNA=subset(scRNA,cells=which(scRNA$seurat_clusters != 25)) #红细胞
scRNA = RenameIdents(scRNA,
                     "0"="T cells", 
                     "1"="Myeloid cells",
                     "2"="T cells",
                     "3"="NK cells",
                     "4"="B cells",
                     "5"="T cells",
                     "6"="T cells",
                     "7"="T cells",
                     "8"="Platelets",
                     "9"="Myeloid cells",
                     "10"="B cells",
                     "11"="Myeloid cells",
                     "12"="Plasma cells",
                     "14"="T cells",
                     "16"="T cells",
                     "17"="Myeloid cells",
                     "18"="Myeloid cells",
                     "19"="T cells",
                     "20"="Myeloid cells",
                     "21"="Myeloid cells",
                     "22"="Myeloid cells",
                     "23"="HSPCs",
                     "24"="B cells",
                     "26"="T cells"
                     )
scRNA$cellType=Idents(scRNA)
saveRDS(scRNA, file="adult.0.5.celltype.Rds")
p <- DimPlot(scRNA, group.by = "cellType", reduction = "umap",
             raster = TRUE,        # 关键：散点栅格化，避免矢量点过多
             pt.size = 1)          # 调整点的大小，看起来更清爽
ggsave("1.1.adult.umap.pdf", plot = p,
       width = 8, height = 6, units = "in",
       dpi = 300, device = "pdf")
markers <- c(
  "CD3D", "CD3E", "CD3G",
  "CD14", "S100A8", "S100A9",
  "NKG7","NCR1","GNLY",
  "CD79A", "CD79B", "MS4A1",
  "PF4", "PPBP", "GP9",
  "JCHAIN", "IGHG1", "SDC1",
  "GATA2", "ETV6", "TYMS"
)
p <- DotPlot(scRNA, features=markers, group.by= "cellType",cols = c("lightgrey", "red"))+
     #旋转横轴标签angle = 90表示将文本旋转90度
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('1.2.adult.marker.pdf', p, width = 10, height = 4)
############### 亚群分析 ###############
cell_type = "Myeloid cells"
aim_cell = subset(scRNA,idents=cell_type)
aim_cell <- NormalizeData(aim_cell, normalization.method = 'LogNormalize', scale.factor = 10000)
aim_cell <- FindVariableFeatures(aim_cell, selection.method = "vst", nfeatures = 2000)
aim_cell <- ScaleData(aim_cell, vars.to.regress = "percent.mt", features = VariableFeatures(aim_cell))
aim_cell <- RunPCA(aim_cell)
pca_embedding <- Embeddings(aim_cell, reduction = "pca")
## 运行Harmony而不依赖Seurat对象
harmony_embedding <- HarmonyMatrix(
  data_mat = pca_embedding,
  meta_data = aim_cell@meta.data,
  vars_use = "orig.ident",
  do_pca = FALSE  # 因为我们已经提供了PCA嵌入
)

## 将结果添加回Seurat对象
aim_cell[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_embedding,
  key = "harmony_",
  assay = "RNA"
)
ggsave('subElbowPlot.pdf', ElbowPlot(aim_cell,ndims = 50), width = 12, height = 10)
#这里是数据选择的最优主成分
pc.num=1:30
#降维聚类
aim_cell <- RunUMAP(aim_cell,reduction="harmony", dims=pc.num,seed.use=123456L)
aim_cell <- RunTSNE(aim_cell,reduction="harmony", dims=pc.num,seed.use=123456L)
aim_cell <- FindNeighbors(aim_cell,reduction="harmony",dims = pc.num,k.param = 20)
aim_cell <- FindClusters(aim_cell,resolution = seq(0.1,0.5,by=0.1))
ggsave('subclustree.pdf', clustree(aim_cell), width = 12, height = 10)
aim_cell$seurat_clusters=aim_cell$RNA_snn_res.0.5
Idents(aim_cell)=aim_cell$seurat_clusters
saveRDS(aim_cell, "sub.adult.M.Rds")
allMarkers <- FindAllMarkers(aim_cell, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(allMarkers, file="allMarkers.aim_cell.M.csv", row.names =F,quote=F)
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 3)) #死细胞
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 6)) #混合细胞
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 8)) #混合细胞
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 12)) #混合细胞
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 14)) #混合细胞
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 15)) #双细胞
aim_cell = RenameIdents(aim_cell,
                     "0"="Classical Monocytes", 
                     "1"="Inflammatory Monocytes",
                     "2"="Non-classical Monocytes",
                     "4"="IFN-responsive Monocytes",
                     "5"="Activated Monocytes",
                     "7"="Tissue-resident Macrophages",
                     "9"="cDC2",
                     "10"="Proliferating Neutrophil Progenitors",
                     "11"="Neutrophils",
                     "13"="Stressed Monocytes"
                     )
aim_cell$subcellType=Idents(aim_cell)
saveRDS(aim_cell, file="sub.adult.M.celltype.Rds")
p=DimPlot(aim_cell,group.by = "subcellType", reduction = "umap", label = FALSE, raster = TRUE, pt.size = 1, cols = subset_colors) #color和heatmap一致
ggsave("2.1.adult.M.umap.pdf", plot = p,
       width = 8, height = 4, units = "in",
       dpi = 300, device = "pdf")

##### T细胞
cell_type = "T cells"
aim_cell = subset(scRNA,idents=cell_type)
aim_cell <- NormalizeData(aim_cell, normalization.method = 'LogNormalize', scale.factor = 10000)
aim_cell <- FindVariableFeatures(aim_cell, selection.method = "vst", nfeatures = 2000)
aim_cell <- ScaleData(aim_cell, vars.to.regress = "percent.mt", features = VariableFeatures(aim_cell))
aim_cell <- RunPCA(aim_cell)
pca_embedding <- Embeddings(aim_cell, reduction = "pca")
## 运行Harmony而不依赖Seurat对象
harmony_embedding <- HarmonyMatrix(
  data_mat = pca_embedding,
  meta_data = aim_cell@meta.data,
  vars_use = "orig.ident",
  do_pca = FALSE  # 因为我们已经提供了PCA嵌入
)

## 将结果添加回Seurat对象
aim_cell[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_embedding,
  key = "harmony_",
  assay = "RNA"
)
#ggsave('subElbowPlot.pdf', ElbowPlot(aim_cell,ndims = 50), width = 12, height = 10)
#这里是数据选择的最优主成分
pc.num=1:30
#降维聚类
aim_cell <- RunUMAP(aim_cell,reduction="harmony", dims=pc.num,seed.use=123456L)
aim_cell <- RunTSNE(aim_cell,reduction="harmony", dims=pc.num,seed.use=123456L)
aim_cell <- FindNeighbors(aim_cell,reduction="harmony",dims = pc.num,k.param = 20)
aim_cell <- FindClusters(aim_cell,resolution = seq(0.1,0.5,by=0.1))
#ggsave('subclustree.pdf', clustree(aim_cell), width = 12, height = 10)
aim_cell$seurat_clusters=aim_cell$RNA_snn_res.0.5
Idents(aim_cell)=aim_cell$seurat_clusters
saveRDS(aim_cell, "sub.adult.T.Rds")
allMarkers <- FindAllMarkers(aim_cell, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(allMarkers, file="allMarkers.aim_cell.T.csv", row.names =F,quote=F)
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 4)) #低质量细胞
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 7)) #低质量细胞
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 11)) #低质量细胞
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 13)) #低质量细胞
aim_cell = RenameIdents(aim_cell,
                     "0"="CD4 Central Memory T cells", 
                     "1"="CD4 Early Activated T cells",
                     "2"="CD8 TEMRA T cells",
                     "3"="CD4 Naive T cells",
                     "5"="CD8 Central Memory T cells",
                     "6"="CD8 TEMRA T cells",
                     "8"="Tregs",
                     "9"="NK-like T cells",
                     "10"="Gamma-delta T cells Vg9Vd2",
                     "12"="MAIT cells",
                     "14"="ISG-high T cells",
                     "15"="CD8 CX3CR1 Effector T cells"
                     )
aim_cell$subcellType=Idents(aim_cell)
saveRDS(aim_cell, file="sub.adult.T.celltype.Rds")
p=DimPlot(aim_cell,group.by = "subcellType", reduction = "umap", label = FALSE, raster = TRUE, pt.size = 1, cols = subset_colors) #color和heatmap一致
ggsave("3.1.adult.T.umap.pdf", plot = p,
       width = 8, height = 6, units = "in",
       dpi = 300, device = "pdf")
###儿童注释
scRN=readRDS("pediatric.0.5.Rds")
scRNA=subset(scRNA,cells=which(scRNA$seurat_clusters != 12)) #红细胞
scRNA = RenameIdents(scRNA,
                     "0"="T cells", 
                     "1"="B cells",
                     "2"="B cells",
                     "3"="T cells",
                     "4"="Myeloid cells",
                     "5"="T cells",
                     "6"="T cells",
                     "7"="T cells",
                     "8"="B cells",
                     "9"="Myeloid cells",
                     "10"="T cells",
                     "11"="Platelets",
                     "13"="Plasma cells",
                     "14"="HSPCs",
                     "15"="Myeloid cells",
                     "16"="T cells",
                     "17"="T cells",
                     "18"="pDCs",
                     "19"="Myeloid cells",
                     "20"="Myeloid cells"
                     )
scRNA$cellType=Idents(scRNA)
saveRDS(scRNA, file="pediatric.0.5.celltype.Rds")
p <- DimPlot(scRNA, group.by = "cellType", reduction = "umap",
             raster = TRUE,        # 关键：散点栅格化，避免矢量点过多
             pt.size = 1)          # 调整点的大小，看起来更清爽
ggsave("1.1.pediatric.umap.pdf", plot = p,
       width = 8, height = 6, units = "in",
       dpi = 300, device = "pdf")
markers <- c(
  "CD79A", "CD79B", "MS4A1",
  "TYMS", "TOP2A", "HMGB2",
  "CD14", "S100A8", "S100A9",
  "NKG7","NCR1","GNLY",
  "LILRA4", "CLEC4C", "IL3RA",
  "JCHAIN", "IGHG1", "SDC1",
  "PF4", "PPBP", "GP9",
  "CD3D", "CD3E", "CD3G"
)
p <- DotPlot(scRNA, features=markers, group.by= "cellType",cols = c("lightgrey", "red"))+
     #旋转横轴标签angle = 90表示将文本旋转90度
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('1.2.pediatric.marker.pdf', p, width = 10, height = 4)

############### 亚群分析 ###############
cell_type = "Myeloid cells"
aim_cell = subset(scRNA,idents=cell_type)
aim_cell <- NormalizeData(aim_cell, normalization.method = 'LogNormalize', scale.factor = 10000)
aim_cell <- FindVariableFeatures(aim_cell, selection.method = "vst", nfeatures = 2000)
aim_cell <- ScaleData(aim_cell, vars.to.regress = "percent.mt", features = VariableFeatures(aim_cell))
aim_cell <- RunPCA(aim_cell)
pca_embedding <- Embeddings(aim_cell, reduction = "pca")
## 运行Harmony而不依赖Seurat对象
harmony_embedding <- HarmonyMatrix(
  data_mat = pca_embedding,
  meta_data = aim_cell@meta.data,
  vars_use = "orig.ident",
  do_pca = FALSE  # 因为我们已经提供了PCA嵌入
)

## 将结果添加回Seurat对象
aim_cell[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_embedding,
  key = "harmony_",
  assay = "RNA"
)
ggsave('subElbowPlot.pdf', ElbowPlot(aim_cell,ndims = 50), width = 12, height = 10)
#这里是数据选择的最优主成分
pc.num=1:30
#降维聚类
aim_cell <- RunUMAP(aim_cell,reduction="harmony", dims=pc.num,seed.use=123456L)
aim_cell <- RunTSNE(aim_cell,reduction="harmony", dims=pc.num,seed.use=123456L)
aim_cell <- FindNeighbors(aim_cell,reduction="harmony",dims = pc.num,k.param = 20)
aim_cell <- FindClusters(aim_cell,resolution = seq(0.1,0.5,by=0.1))
ggsave('subclustree.pdf', clustree(aim_cell), width = 12, height = 10)
aim_cell$seurat_clusters=aim_cell$RNA_snn_res.0.5
Idents(aim_cell)=aim_cell$seurat_clusters
saveRDS(aim_cell, "sub.pediatric.M.Rds")
allMarkers <- FindAllMarkers(aim_cell, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(allMarkers, file="allMarkers.aim_cell.M.csv", row.names =F,quote=F)
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 4)) #死细胞
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 6)) #双细胞
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 8)) #双细胞
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 9)) #死细胞
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 11)) #红细胞
aim_cell = RenameIdents(aim_cell,
                     "0"="Tissue-resident Macrophages", 
                     "1"="Neutrophils",
                     "2"="Inflammatory Monocytes",
                     "3"="IFN-responsive Monocytes",
                     "5"="Proliferating Neutrophil Progenitors",
                     "7"="Stressed Monocytes",
                     "10"="cDC2"
                     )
aim_cell$subcellType=Idents(aim_cell)
saveRDS(aim_cell, file="sub.pediatric.M.celltype.Rds")
p=DimPlot(aim_cell,group.by = "subcellType", reduction = "umap", label = FALSE, raster = TRUE, pt.size = 1.5, cols = subset_colors) #color和heatmap一致
ggsave("2.1.pediatric.M.umap.pdf", plot = p,
       width = 8, height = 5, units = "in",
       dpi = 300, device = "pdf")
##### T细胞
cell_type = "T cells"
aim_cell = subset(scRNA,idents=cell_type)
aim_cell <- NormalizeData(aim_cell, normalization.method = 'LogNormalize', scale.factor = 10000)
aim_cell <- FindVariableFeatures(aim_cell, selection.method = "vst", nfeatures = 2000)
aim_cell <- ScaleData(aim_cell, vars.to.regress = "percent.mt", features = VariableFeatures(aim_cell))
aim_cell <- RunPCA(aim_cell)
pca_embedding <- Embeddings(aim_cell, reduction = "pca")
## 运行Harmony而不依赖Seurat对象
harmony_embedding <- HarmonyMatrix(
  data_mat = pca_embedding,
  meta_data = aim_cell@meta.data,
  vars_use = "orig.ident",
  do_pca = FALSE  # 因为我们已经提供了PCA嵌入
)

## 将结果添加回Seurat对象
aim_cell[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_embedding,
  key = "harmony_",
  assay = "RNA"
)
ggsave('subElbowPlot.pdf', ElbowPlot(aim_cell,ndims = 50), width = 12, height = 10)
#这里是数据选择的最优主成分
pc.num=1:30
#降维聚类
aim_cell <- RunUMAP(aim_cell,reduction="harmony", dims=pc.num,seed.use=123456L)
aim_cell <- RunTSNE(aim_cell,reduction="harmony", dims=pc.num,seed.use=123456L)
aim_cell <- FindNeighbors(aim_cell,reduction="harmony",dims = pc.num,k.param = 20)
aim_cell <- FindClusters(aim_cell,resolution = seq(0.1,0.5,by=0.1))
ggsave('subclustree.pdf', clustree(aim_cell), width = 12, height = 10)
aim_cell$seurat_clusters=aim_cell$RNA_snn_res.0.5
Idents(aim_cell)=aim_cell$seurat_clusters
saveRDS(aim_cell, "sub.pediatric.T.Rds")
allMarkers <- FindAllMarkers(aim_cell, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(allMarkers, file="allMarkers.aim_cell.T.csv", row.names =F,quote=F)

# 将NK细胞还到细胞大类上
nk_barcodes <- colnames(aim_cell)[aim_cell$seurat_clusters == 6]
# 更新 scRNA 中的 cellType 列
scRNA$cellType <- as.character(scRNA$cellType)
scRNA$cellType[colnames(scRNA) %in% nk_barcodes] <- "NK cells"
scRNA$cellType <- factor(scRNA$cellType)  # 可选
# 检查更新结果
table(scRNA$cellType)

aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 4)) #低质量细胞
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 10)) #血小板
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 11)) #B细胞和单核细胞混合细胞
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 13)) #浆细胞
aim_cell=subset(aim_cell,cells=which(aim_cell$seurat_clusters != 14)) #HSPCs
aim_cell = RenameIdents(aim_cell,
                     "0"="CD4 Naive T cells", 
                     "1"="CD8 Central Memory T cells",
                     "2"="CD4 Central Memory T cells",
                     "3"="CD4 Early Activated T cells",
                     "5"="CD8 TEMRA T cells",
                     "6"="CD8 TEMRA T cells",
                     "7"="Gamma-delta T cells Vg9Vd2",
                     "8"="Gamma-delta T cells Vd1",
                     "9"="Tregs",
                     "12"="ISG-high T cells"
                     )
aim_cell$subcellType=Idents(aim_cell)
saveRDS(aim_cell, file="sub.pediatric.T.celltype.Rds")
p=DimPlot(aim_cell,group.by = "subcellType", reduction = "umap", label = FALSE, raster = TRUE, pt.size = 1.5, cols = subset_colors) #color和heatmap一致
ggsave("3.1.pediatric.T.umap.pdf", plot = p,
       width = 8, height = 5, units = "in",
       dpi = 600, device = "pdf")




