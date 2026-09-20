#拟时序分析
library(slingshot)
library(scuttle)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(fields) 
library(ggrastr)
#儿童M细胞 
sce <- as.SingleCellExperiment(aim_cell)
sce <- slingshot(sce, 
                 clusterLabels = "subcellType",  # 你的注释列名
                 reducedDim = "UMAP",         # 使用UMAP降维
                 start.clus = "Inflammatory Monocytes",          # 指定起始细胞类型（根据你的生物学知识）
                 approx_points = 100)         # 轨迹线平滑度
print(slingLineages(sce))
saveRDS(sce, file="pediatric.M.slingshot.Rds")

# ---------- 1. 提取 PCA 坐标 ----------
pca <- reducedDims(sce)$PCA
pca_df <- data.frame(Dim1 = pca[, 1], Dim2 = pca[, 2])
pca_df$cell <- colnames(sce)
# ---------- 2. 提取 Lineage2 伪时间 ----------
pt <- sce$slingPseudotime_2          # 如列名不同请相应修改
pca_df$pseudotime <- pt
# ---------- 3. 提取 Lineage2 曲线坐标 ----------
crv3 <- slingCurves(sce)[["Lineage2"]]
curve_df <- data.frame(Dim1 = crv3$s[, 1], Dim2 = crv3$s[, 2])
# ---------- 4. 渐变色（蓝→红，早→晚） ----------
colors <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(100)
plotcol <- colors[cut(pt, breaks = 100, labels = FALSE)]

pdf("2.7.M.slingshot.pdf", width = 7, height = 6, bg = "white")
# 主图：PCA 散点 + 轨迹曲线
plot(pca_df$Dim1, pca_df$Dim2,
     col = plotcol, pch = 16, cex = 0.6, asp = 1,
     xlab = "PC1", ylab = "PC2",
     main = "",
     bty = "l", mgp = c(2, 0.5, 0))

# 叠加曲线（黑色加粗）
lines(curve_df$Dim1, curve_df$Dim2, lwd = 3, col = "black")
# 伪时间图例（右侧色条）
image.plot(legend.only = TRUE,
           zlim = range(pt, na.rm = TRUE),
           col = colors,
           legend.lab = "Pseudotime",
           legend.line = 2.5,
           legend.mar = 6)
dev.off()

###########多条轨迹
curves <- slingCurves(sce)
curve_df <- map_dfr(seq_along(curves), function(i) {
  c <- curves[[i]]
  data.frame(
    UMAP1 = c$s[c$ord, 1],
    UMAP2 = c$s[c$ord, 2],
    lineage = i
  )
})

# 2. 提取细胞UMAP坐标和伪时间矩阵，确定每个细胞的主谱系
umap_df <- as.data.frame(reducedDim(sce, "UMAP")) %>%
  setNames(c("UMAP1", "UMAP2"))

pseudo_mat <- slingPseudotime(sce)
branch <- apply(pseudo_mat, 1, function(x) {
  non_na <- which(!is.na(x))
  if (length(non_na) == 1) non_na else 0L
})
umap_df$lineage <- branch

# 3. 为每个谱系分配颜色渐变（viridis调色板循环）
palettes <- c("plasma", "magma", "inferno", "cividis")
umap_df$color <- NA
for (i in seq_len(ncol(pseudo_mat))) {
  idx <- umap_df$lineage == i
  if (!any(idx)) next
  pt_norm <- (pseudo_mat[idx, i] - min(pseudo_mat[, i], na.rm = TRUE)) /
             diff(range(pseudo_mat[, i], na.rm = TRUE))
  pal <- viridis::viridis_pal(option = palettes[(i-1) %% length(palettes) + 1])(100)
  umap_df$color[idx] <- pal[pmax(1, pmin(100, round(pt_norm * 99) + 1))]
}
umap_df$color[umap_df$lineage == 0] <- "#0d0887"  # 根细胞颜色

# 4. 绘图
p=ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(aes(color = color), size = 0.5, alpha = 0.9) +
  scale_color_identity() +
  geom_path(data = curve_df, aes(group = lineage), color = "black", linewidth = 1.2) +
  theme_classic()
ggsave('2.7.M.slingshot.pdf', p, width = 8, height = 6)

# cellchat细胞通讯
library(CellChat)
cellchat <- createCellChat(object = aim_cell,
                           meta = aim_cell@meta.data,
                           group.by = "subcellType")
CellChatDB <- CellChatDB.human  
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "2.7.M.cell_communications.csv")
#计算每个信号通路相关的所有配体-受体相互作用的通信结果
cellchat <- computeCommunProbPathway(cellchat)
#计算整合的细胞类型之间通信结果
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

p=netVisual_bubble(cellchat, sources.use = c(3),targets.use = c(1,2,4,5,6,7), remove.isolate = FALSE)
ggsave("2.7.M.Cellchat.pdf",p,width = 5, height = 8)

#儿童T细胞
sce <- as.SingleCellExperiment(aim_cell)
sce <- slingshot(sce, 
                 clusterLabels = "subcellType",  # 你的注释列名
                 reducedDim = "UMAP",         # 使用UMAP降维
                 start.clus = "CD4 Naive T cells",                   
                 approx_points = 100)         # 轨迹线平滑度
print(slingLineages(sce))
saveRDS(sce, file="pediatric.T.slingshot.Rds")
curves <- slingCurves(sce)
curve_df <- map_dfr(seq_along(curves), function(i) {
  c <- curves[[i]]
  data.frame(
    UMAP1 = c$s[c$ord, 1],
    UMAP2 = c$s[c$ord, 2],
    lineage = i
  )
})

# 2. 提取细胞UMAP坐标和伪时间矩阵，确定每个细胞的主谱系
umap_df <- as.data.frame(reducedDim(sce, "UMAP")) %>%
  setNames(c("UMAP1", "UMAP2"))

pseudo_mat <- slingPseudotime(sce)
branch <- apply(pseudo_mat, 1, function(x) {
  non_na <- which(!is.na(x))
  if (length(non_na) == 1) non_na else 0L
})
umap_df$lineage <- branch

# 3. 为每个谱系分配颜色渐变（viridis调色板循环）
palettes <- c("plasma", "magma", "inferno", "cividis")
umap_df$color <- NA
for (i in seq_len(ncol(pseudo_mat))) {
  idx <- umap_df$lineage == i
  if (!any(idx)) next
  pt_norm <- (pseudo_mat[idx, i] - min(pseudo_mat[, i], na.rm = TRUE)) /
             diff(range(pseudo_mat[, i], na.rm = TRUE))
  pal <- viridis::viridis_pal(option = palettes[(i-1) %% length(palettes) + 1])(100)
  umap_df$color[idx] <- pal[pmax(1, pmin(100, round(pt_norm * 99) + 1))]
}
umap_df$color[umap_df$lineage == 0] <- "#0d0887"  # 根细胞颜色

# 4. 绘图
p=ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(aes(color = color), size = 0.5, alpha = 0.9) +
  scale_color_identity() +
  geom_path(data = curve_df, aes(group = lineage), color = "black", linewidth = 1.2) +
  theme_classic()
ggsave('3.7.T.slingshot.pdf', p, width = 8, height = 6)
# cellchat细胞通讯
library(CellChat)
cellchat <- createCellChat(object = aim_cell,
                           meta = aim_cell@meta.data,
                           group.by = "subcellType")
CellChatDB <- CellChatDB.human  
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "3.7.T.cell_communications.csv")
#计算每个信号通路相关的所有配体-受体相互作用的通信结果
cellchat <- computeCommunProbPathway(cellchat)
#计算整合的细胞类型之间通信结果
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

p=netVisual_bubble(cellchat, sources.use = c(2),targets.use = c(1,3,4,5,6,7,8,9), remove.isolate = FALSE)
ggsave("3.7.T.Cellchat.pdf",p,width = 5, height = 4)


#成人M细胞
sce <- as.SingleCellExperiment(aim_cell)
sce <- slingshot(sce, 
                 clusterLabels = "subcellType",  # 你的注释列名
                 reducedDim = "UMAP",         # 使用UMAP降维
                 start.clus = "cDC2",          # 指定起始细胞类型（根据你的生物学知识）
                 approx_points = 100)         # 轨迹线平滑度
print(slingLineages(sce))
saveRDS(sce, file="adult.M.slingshot.Rds")
# ---------- 1. 提取 PCA 坐标 ----------
pca <- reducedDims(sce)$PCA
pca_df <- data.frame(Dim1 = pca[, 1], Dim2 = pca[, 2])
pca_df$cell <- colnames(sce)
# ---------- 2. 提取 Lineage3 伪时间 ----------
pt <- sce$slingPseudotime_3          # 如列名不同请相应修改
pca_df$pseudotime <- pt
# ---------- 3. 提取 Lineage3 曲线坐标 ----------
crv3 <- slingCurves(sce)[["Lineage3"]]
curve_df <- data.frame(Dim1 = crv3$s[, 1], Dim2 = crv3$s[, 2])
# ---------- 4. 渐变色（蓝→红，早→晚） ----------
colors <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(100)
plotcol <- colors[cut(pt, breaks = 100, labels = FALSE)]

pdf("2.7.M.slingshot.pdf", width = 7, height = 6, bg = "white")
# 主图：PCA 散点 + 轨迹曲线
plot(pca_df$Dim1, pca_df$Dim2,
     col = plotcol, pch = 16, cex = 0.6, asp = 1,
     xlab = "PC1", ylab = "PC2",
     main = "",
     bty = "l", mgp = c(2, 0.5, 0))

# 叠加曲线（黑色加粗）
lines(curve_df$Dim1, curve_df$Dim2, lwd = 3, col = "black")
# 伪时间图例（右侧色条）
image.plot(legend.only = TRUE,
           zlim = range(pt, na.rm = TRUE),
           col = colors,
           legend.lab = "Pseudotime",
           legend.line = 2.5,
           legend.mar = 6)
dev.off()


# cellchat细胞通讯
library(CellChat)
cellchat <- createCellChat(object = aim_cell,
                           meta = aim_cell@meta.data,
                           group.by = "subcellType")
CellChatDB <- CellChatDB.human  
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "2.7.M.cell_communications.csv")
#计算每个信号通路相关的所有配体-受体相互作用的通信结果
cellchat <- computeCommunProbPathway(cellchat)
#计算整合的细胞类型之间通信结果
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

p=netVisual_bubble(cellchat, sources.use = c(7),targets.use = c(1,2,4,5,6,8,9,10), remove.isolate = FALSE)
ggsave("2.7.M.Cellchat.pdf",p,width = 5, height = 8)

#成人T细胞
sce <- as.SingleCellExperiment(aim_cell)
sce <- slingshot(sce, 
                 clusterLabels = "subcellType",  # 你的注释列名
                 reducedDim = "UMAP",         # 使用UMAP降维
                 start.clus = "CD4 Naive T cells",          # 指定起始细胞类型（根据你的生物学知识）
                 approx_points = 100)         # 轨迹线平滑度
print(slingLineages(sce))
saveRDS(sce, file="adult.T.slingshot.Rds")
# 1. 提取 PCA 坐标
lin1 <- slingLineages(sce)[["Lineage1"]]
sel  <- sce$subcellType %in% lin1                    
pt_all <- slingPseudotime(sce, na = FALSE)[,         
            match("Lineage1", colnames(slingPseudotime(sce)))]
names(pt_all) <- colnames(sce)                       

pca <- reducedDims(sce)$PCA
pca_df <- data.frame(Dim1 = pca[sel, 1], Dim2 = pca[sel, 2])      
pca_df$cell       <- colnames(sce)[sel]                            
pca_df$celltype   <- droplevels(factor(sce$subcellType[sel]))      
pca_df$pseudotime <- pt_all[sel]                                   

centroids <- aggregate(cbind(Dim1, Dim2) ~ celltype, data = pca_df, FUN = mean)

p <- ggplot(pca_df, aes(Dim1, Dim2, color = pseudotime)) +
  rasterise(geom_point(size = 0.4), dpi = 300, scale = 1) +
  scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
  geom_text(data = centroids,
            aes(x = Dim1, y = Dim2, label = celltype),
            inherit.aes = FALSE,
            size = 3, fontface = "bold") +
  coord_fixed() +
  theme_classic()
ggsave("3.7.T.slingshot.pdf", p, width = 8, height = 7, dpi = 300, device = "pdf")

# cellchat细胞通讯
library(CellChat)
cellchat <- createCellChat(object = aim_cell,
                           meta = aim_cell@meta.data,
                           group.by = "subcellType")
CellChatDB <- CellChatDB.human  
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "3.7.T.cell_communications.csv")
#计算每个信号通路相关的所有配体-受体相互作用的通信结果
cellchat <- computeCommunProbPathway(cellchat)
#计算整合的细胞类型之间通信结果
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

p=netVisual_bubble(cellchat, sources.use = c(3),targets.use = c(1,2,4,5,6,7,8,9,10,11), remove.isolate = FALSE)
ggsave("3.7.T.Cellchat.pdf",p,width = 5, height = 4)