###### 结合临床信息分析
aim_cell=readRDS("sub.adult.M.celltype.Rds")
colnames(aim_cell@meta.data)[colnames(aim_cell@meta.data) == "Gram.egative.bacteria"] <- "Gram.negative.bacteria"
library(tidyverse)
library(broom)
library(ggplot2)
### 创建临床表型表格
clinical <- aim_cell@meta.data %>%
    dplyr::distinct(orig.ident, .keep_all = TRUE) %>%
    dplyr::select(
        donor = orig.ident,
        Age,
        Gender,
        SOFA = pSOFA.score,
        death_28 = X28day,
        Gram_pos = Gram.positive.bacteria,
        Gram_neg = Gram.negative.bacteria,  
        Virus,
        Fungi,
        Mixed = Mixed.infection,
        Unknown
    ) %>%
    mutate(
        inf_group = case_when(
            Gram_pos == 1 & Gram_neg == 0 & Virus == 0 & Fungi == 0 & Mixed == 0 ~ "Gram_pos",
            Gram_neg == 1 & Gram_pos == 0 & Virus == 0 & Fungi == 0 & Mixed == 0 ~ "Gram_neg",
            Virus == 1 & Gram_pos == 0 & Gram_neg == 0 & Fungi == 0 & Mixed == 0 ~ "Virus",
            Fungi == 1 & Gram_pos == 0 & Gram_neg == 0 & Virus == 0 & Mixed == 0 ~ "Fungi",
            Mixed == 1 ~ "Mixed",
            Unknown == 1 ~ "Unknown",
            TRUE ~ "Other"
        )
    ) %>%
    dplyr::select(-starts_with("Gram_"), -Virus, -Fungi, -Mixed, -Unknown)
##adult
clinical <- aim_cell@meta.data %>%
    dplyr::distinct(orig.ident, .keep_all = TRUE) %>%
    dplyr::select(
        donor = orig.ident,
        Age,
        Gender,
        SOFA = SOFA.score,
        death_28 = X28day,
        Gram_pos = Gram.positive.bacteria,
        Gram_neg = Gram.negative.bacteria,  
        Virus,
        Fungi,
        Unknown
    ) %>%
    mutate(
        inf_group = case_when(
            Gram_pos == 1 & Gram_neg == 0 & Virus == 0 & Fungi == 0 ~ "Gram_pos",
            Gram_neg == 1 & Gram_pos == 0 & Virus == 0 & Fungi == 0 ~ "Gram_neg",
            Virus == 1 & Gram_pos == 0 & Gram_neg == 0 & Fungi == 0 ~ "Virus",
            Fungi == 1 & Gram_pos == 0 & Gram_neg == 0 & Virus == 0 ~ "Fungi",
            Unknown == 1 ~ "Unknown",
            TRUE ~ "Other"
        )
    ) %>%
    dplyr::select(-starts_with("Gram_"), -Virus, -Fungi, -Unknown)
### 比例关联
prop_tab <- table(aim_cell$orig.ident, aim_cell$subcellType)
prop_mat <- prop.table(prop_tab, margin = 1)
prop_df <- as.data.frame.matrix(prop_mat)
prop_df$donor <- rownames(prop_df)
dat <- left_join(prop_df, clinical, by = "donor")
# 亚群名
subclusters <- setdiff(colnames(prop_mat), "donor")

# SOFA评分：逐亚群自动离群值检测 + 相关性分析
source("3.2.Clean&Cor&Plot.R")
# ---------- 对所有亚群执行自动清洗和相关性分析 ----------
results_all <- map_dfr(subclusters, function(sc) {
  auto_clean_and_cor(dat, sc, clinical_var = "SOFA")
})

# 计算清洗后的 FDR
results_all$fdr_clean <- p.adjust(results_all$p_clean, method = "BH")

# 按清洗后 p 值排序
results_all <- results_all %>% arrange(p_clean)

# 打印完整结果
write.csv(results_all, file="2.3.M.SOFA.csv", row.names =F,quote=F)
print("=== 自动离群值清洗后的相关性结果 ===")
print(results_all, n = Inf)

# ---------- 3. 相关性分析结果可视化 ----------
# 提取两个目标亚群的清洗后数据
target_clusters <- "cDC2"

cleaned_list <- lapply(target_clusters, function(sc) {
  auto_clean_and_plot(dat, sc, clinical_var = "SOFA")
})

# 合并为一个数据框
plot_data <- do.call(rbind, cleaned_list)
colnames(plot_data) <- c("Proportion", "SOFA", "Subcluster")
# ============================================================
# 3. 从 results_all 中提取清洗后的相关系数
# ============================================================
cor_labels <- results_all %>%
  filter(subcluster %in% target_clusters) %>%
  mutate(
    Subcluster = subcluster,  # 重命名为与分面变量一致
    label = paste0("rho = ", round(rho_clean, 3), 
                   ", p = ", format(p_clean, digits = 3))
  ) %>%
  select(Subcluster, label)
# 绘图
p <- ggplot(plot_data, aes(x = SOFA, y = Proportion)) +
  geom_point(size = 3, alpha = 0.8, color = "#2c3e50") +
  geom_smooth(method = "lm", se = TRUE, color = "#e74c3c", fill = "#e74c3c", alpha = 0.2) +
  facet_wrap(~ Subcluster, scales = "free_y", ncol = 2) +
  geom_text(data = cor_labels, 
            aes(label = label),
            x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
            size = 4, color = "black", fontface = "italic") + 
  labs(
    x = "SOFA Score",
    y = "Proportion of cells"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )
ggsave("2.3.M.SOFA.pdf", p, width = 7, height = 5)

# 各个亚群之间的相关性
library(pheatmap)
library(RColorBrewer)
library(Hmisc)
# 构建相关性矩阵
res <- rcorr(as.matrix(dat[, subclusters]), type = "spearman")
cor_matrix <- res$r   # 相关系数矩阵
p_mat <- res$P        # P 值矩阵
# 构造带换行的文本矩阵（如上所示）
text_matrix <- matrix("", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    r <- cor_matrix[i, j]
    p <- p_mat[i, j]
    star <- ""
    if (!is.na(p)) {                    
      if (p < 0.001) star <- "***"
      else if (p < 0.01) star <- "**"
      else if (p < 0.05) star <- "*"
    }
    if (star != "") {
      text_matrix[i, j] <- sprintf("%.2f\n%s", r, star)
    } else {
      text_matrix[i, j] <- sprintf("%.2f\n ", r)
    }
  }
}
rownames(text_matrix) <- rownames(cor_matrix)
colnames(text_matrix) <- colnames(cor_matrix)
# 绘制热图
pheatmap(cor_matrix,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         display_numbers = text_matrix,
         number_color = "black",
         fontsize_number = 10,
         fontsize_row = 10,
         fontsize_col = 10,
         main = "Correlation Heatmap between T Subclusters",
         angle_col = 315,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         legend = FALSE,
         filename = "3.4.T.cor_subcellType.pdf",
         width = 8,
         height = 6)