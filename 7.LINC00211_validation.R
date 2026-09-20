suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
})

# --- 路径配置 ---
ped_expr_path <- "pediatric/2.6.bulk_tmm_log2cpm.csv"
ped_clin_path <- "pediatric/2.6.bulk_clinical.csv"
adult_expr_path <- "adult/2.6.bulk_tmm_log2cpm.csv"
adult_clin_path <- "adult/2.6.bulk_clinical.csv"

out_dir <- "LINC00211_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 加载数据 ---
cat("=== Loading data ===\n")
ped_expr <- read.csv(ped_expr_path, row.names = 1, check.names = FALSE)
ped_clin <- read.csv(ped_clin_path, stringsAsFactors = FALSE)
adult_expr <- read.csv(adult_expr_path, row.names = 1, check.names = FALSE)
adult_clin <- read.csv(adult_clin_path, stringsAsFactors = FALSE)

cat("Pediatric:", nrow(ped_expr), "genes x", ncol(ped_expr), "samples |",
    nrow(ped_clin), "clinical records\n")
cat("Adult:", nrow(adult_expr), "genes x", ncol(adult_expr), "samples |",
    nrow(adult_clin), "clinical records\n")

# 确认列对齐
stopifnot(all(colnames(ped_expr) == ped_clin$donor))
stopifnot(all(colnames(adult_expr) == adult_clin$donor))

# 儿童髓系 consistent genes (16)
ped_myeloid_genes <- c("WNK1","CCT7","HDGF","FBXO7","TNIP2","UBE3C","CBX1","TRPV2",
                       "VGLL4","SUN2","MVD","IARS2","U2AF2","DCXR","CSNK2A2","MPHOSPH6")
# 儿童T consistent genes (3)
ped_tcell_genes <- c("UFD1","IMMT","THOP1")
# 成人髓系 consistent genes (3)
adult_myeloid_genes <- c("SOD1","MGA","MTMR9")
# 成人T consistent genes (5)
adult_tcell_genes <- c("UPF3A","SOCS1","CYB561D1","ND4L","PLIN2")
# 跨模块枢纽
il16_gene <- c("IL16")
# 轨迹通讯配体（CellChat 中"目标亚群→选中轨迹细胞"实际参与通讯的配体，与共表达输入一致）
# 儿童髓系 IM→Lineage2：CCL3/CCL3L1/MIF/IL16/RETN/NAMPT/ANXA1/GRN/LGALS9
# 成人髓系 cDC2→Lineage3：MIF/IL16/ANXA1/GAS6/GRN/LGALS9
# 成人T CD8 TEMRA：TNF/MIF/IL16/GZMA
ligand_genes <- c("CCL3","CCL3L1","MIF","RETN","NAMPT","ANXA1","GRN","LGALS9",
                  "GAS6","TNF","GZMA")
pirat_targets <- c("S100A8","S100A9")

myeloid_genes <- c(ped_myeloid_genes, adult_myeloid_genes)   # 髓系（正面验证）
tcell_genes  <- c(ped_tcell_genes, adult_tcell_genes)        # T 细胞（阴性对照）

all_target_genes <- unique(c(myeloid_genes, tcell_genes, il16_gene,
                             pirat_targets, ligand_genes))

# 基因模块分类（用于结果标注）
module_cat <- c(
  WNK1="Ped.Myeloid", CCT7="Ped.Myeloid", HDGF="Ped.Myeloid", FBXO7="Ped.Myeloid",
  TNIP2="Ped.Myeloid", UBE3C="Ped.Myeloid", CBX1="Ped.Myeloid", TRPV2="Ped.Myeloid",
  VGLL4="Ped.Myeloid", SUN2="Ped.Myeloid", MVD="Ped.Myeloid", IARS2="Ped.Myeloid",
  U2AF2="Ped.Myeloid", DCXR="Ped.Myeloid", CSNK2A2="Ped.Myeloid", MPHOSPH6="Ped.Myeloid",
  UFD1="Ped.T", IMMT="Ped.T", THOP1="Ped.T",
  SOD1="Adult.Myeloid", MGA="Adult.Myeloid", MTMR9="Adult.Myeloid",
  UPF3A="Adult.T", SOCS1="Adult.T", CYB561D1="Adult.T", ND4L="Adult.T", PLIN2="Adult.T",
  IL16="IL16",
  S100A8="PIRAT1.Target", S100A9="PIRAT1.Target",
  CCL3="Ligand", CCL3L1="Ligand", MIF="Ligand", RETN="Ligand", NAMPT="Ligand",
  ANXA1="Ligand", GRN="Ligand", LGALS9="Ligand", GAS6="Ligand", TNF="Ligand", GZMA="Ligand"
)
# 分析分组（摘要/着色用）：髓系正面 / T细胞阴性对照 / 机制 / 通讯 / 跨模块
analysis_group <- function(g) {
  if (g %in% c(ped_myeloid_genes, adult_myeloid_genes)) "Myeloid (positive)"
  else if (g %in% c(ped_tcell_genes, adult_tcell_genes)) "T cell (negative control)"
  else if (g %in% pirat_targets) "PIRAT1 Target (mechanism)"
  else if (g %in% ligand_genes) "CellChat Ligand"
  else "IL16 (cross-module)"
}

# 确认基因存在
for (g in c("LINC00211", all_target_genes)) {
  ped_ok <- g %in% rownames(ped_expr)
  adult_ok <- g %in% rownames(adult_expr)
  if (!ped_ok || !adult_ok) {
    cat("WARNING:", g, "- Pediatric:", ped_ok, "Adult:", adult_ok, "\n")
  }
}

# --- LINC00211 表达概况 ---
cat("\n=== LINC00211 Expression Overview ===\n")
for (label in c("Pediatric","Adult")) {
  expr <- if (label == "Pediatric") ped_expr else adult_expr
  clin <- if (label == "Pediatric") ped_clin else adult_clin
  vals <- as.numeric(expr["LINC00211", ])
  cat(sprintf("%s: range [%.2f, %.2f] | median %.2f | IQR [%.2f, %.2f] | >0: %d/%d\n",
              label, min(vals), max(vals), median(vals),
              quantile(vals, 0.25), quantile(vals, 0.75),
              sum(vals > 0), length(vals)))
}

# Part 1: 临床关联分析（SOFA 梯度 + 28 天死亡）
cat("\n=== Part 1: Clinical Association ===\n")

clinical_results <- data.frame()

for (label in c("Pediatric","Adult")) {
  expr <- if (label == "Pediatric") ped_expr else adult_expr
  clin <- if (label == "Pediatric") ped_clin else adult_clin
  n <- nrow(clin)

  lnc_vals <- as.numeric(expr["LINC00211", ])

  # --- SOFA Spearman ---
  sofa_cor <- cor.test(lnc_vals, clin$SOFA, method = "spearman", exact = FALSE)
  cat(sprintf("[%s] LINC00211 ~ SOFA: rho = %.3f, p = %.4f (n=%d)\n",
              label, sofa_cor$estimate, sofa_cor$p.value, n))

  # --- 28-day death ---
  death_col <- "death_28"
  death_vals <- clin[[death_col]]
  death_valid <- !is.na(death_vals)
  if (sum(death_valid) > 0 && length(unique(death_vals[death_valid])) > 1) {
    death_test <- wilcox.test(lnc_vals[death_valid] ~ factor(death_vals[death_valid]))
    cat(sprintf("[%s] LINC00211 ~ 28-day death: W = %.1f, p = %.4f (alive=%d, dead=%d)\n",
                label, death_test$statistic, death_test$p.value,
                sum(death_vals[death_valid] == 0), sum(death_vals[death_valid] == 1)))

    # Logistic regression (control Age + Gender)
    if (sum(death_vals[death_valid] == 1) >= 3) {
      logi_df <- data.frame(
        LINC00211 = lnc_vals[death_valid],
        death = death_vals[death_valid],
        Age = clin$Age[death_valid],
        Gender = factor(clin$Gender[death_valid])
      )
      logi_fit <- glm(death ~ LINC00211 + Age + Gender, data = logi_df, family = binomial())
      logi_sum <- summary(logi_fit)
      or <- exp(coef(logi_fit)["LINC00211"])
      or_ci <- exp(confint(logi_fit)["LINC00211", ])
      cat(sprintf("[%s] Logistic OR = %.3f [%.3f, %.3f], p = %.4f\n",
                  label, or, or_ci[1], or_ci[2], logi_sum$coefficients["LINC00211","Pr(>|z|)"]))
    }
  }

  # 记录结果
  clinical_results <- rbind(clinical_results, data.frame(
    cohort = label,
    n = n,
    sofa_rho = sofa_cor$estimate,
    sofa_p = sofa_cor$p.value,
    death_events = sum(death_vals == 1, na.rm = TRUE),
    stringsAsFactors = FALSE
  ))
}

write.csv(clinical_results, file.path(out_dir, "LINC00211_clinical_association.csv"), row.names = FALSE)

# Part 2: 共表达特异性分析
cat("\n=== Part 2: Co-expression Specificity ===\n")

coexpr_all <- data.frame()

for (label in c("Pediatric","Adult")) {
  expr <- if (label == "Pediatric") ped_expr else adult_expr
  clin <- if (label == "Pediatric") ped_clin else adult_clin

  lnc_vals <- as.numeric(expr["LINC00211", ])

  for (gene in all_target_genes) {
    if (!gene %in% rownames(expr)) next
    gene_vals <- as.numeric(expr[gene, ])

    # Spearman
    sp <- cor.test(lnc_vals, gene_vals, method = "spearman", exact = FALSE)
    # Pearson
    pe <- cor.test(lnc_vals, gene_vals, method = "pearson")

    # 分类
    category <- analysis_group(gene)
    module <- module_cat[[gene]]

    coexpr_all <- rbind(coexpr_all, data.frame(
      cohort = label,
      gene = gene,
      module = module,
      category = category,
      spearman_rho = sp$estimate,
      spearman_p = sp$p.value,
      pearson_r = pe$estimate,
      pearson_p = pe$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

coexpr_all$spearman_fdr <- p.adjust(coexpr_all$spearman_p, method = "BH")
coexpr_all$significant <- coexpr_all$spearman_p < 0.05 & abs(coexpr_all$spearman_rho) > 0.3

write.csv(coexpr_all, file.path(out_dir, "LINC00211_coexpression.csv"), row.names = FALSE)

# Part 3: 偏相关（控制 SOFA）
cat("\n=== Part 3: Partial Correlation (controlling SOFA) ===\n")

partial_results <- data.frame()

for (label in c("Pediatric","Adult")) {
  expr <- if (label == "Pediatric") ped_expr else adult_expr
  clin <- if (label == "Pediatric") ped_clin else adult_clin

  lnc_vals <- as.numeric(expr["LINC00211", ])
  sofa <- clin$SOFA

  for (gene in all_target_genes) {
    if (!gene %in% rownames(expr)) next
    gene_vals <- as.numeric(expr[gene, ])

    # 偏相关: LINC00211 ~ Gene + SOFA, 取 LINC00211 残差 vs Gene 残差
    df_cor <- data.frame(lnc = lnc_vals, gene = gene_vals, sofa = sofa)
    df_cor <- df_cor[complete.cases(df_cor), ]

    if (nrow(df_cor) < 5) next

    # 回归去除 SOFA
    lnc_res <- residuals(lm(lnc ~ sofa, data = df_cor))
    gene_res <- residuals(lm(gene ~ sofa, data = df_cor))

    sp <- cor.test(lnc_res, gene_res, method = "spearman", exact = FALSE)
    pe <- cor.test(lnc_res, gene_res, method = "pearson")

    partial_results <- rbind(partial_results, data.frame(
      cohort = label,
      gene = gene,
      module = module_cat[[gene]],
      category = analysis_group(gene),
      partial_rho = sp$estimate,
      partial_p = sp$p.value,
      partial_r = pe$estimate,
      partial_p_pearson = pe$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

partial_results$partial_fdr <- p.adjust(partial_results$partial_p, method = "BH")

write.csv(partial_results, file.path(out_dir, "LINC00211_partial_cor.csv"), row.names = FALSE)

# 可视化
cat("\n=== Part 5: Visualization ===\n")

# --- Fig 1: SOFA 散点图 ---
plot_list <- list()
for (label in c("Pediatric","Adult")) {
  expr <- if (label == "Pediatric") ped_expr else adult_expr
  clin <- if (label == "Pediatric") ped_clin else adult_clin
  lnc_vals <- as.numeric(expr["LINC00211", ])

  p <- ggplot(data.frame(LINC00211 = lnc_vals, SOFA = clin$SOFA,
                         death = factor(clin$death_28)),
              aes(x = SOFA, y = LINC00211, color = death)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40", linewidth = 0.5) +
    stat_cor(method = "spearman", label.x.npc = 0.05, label.y.npc = 0.95,
             aes(label = ..r.label..), color = "black", size = 3.5) +
    scale_color_manual(values = c("0" = "#4DBBD5", "1" = "#E64B35"),
                       labels = c("0" = "Alive", "1" = "Deceased"),
                       name = "28-day") +
    labs(title = paste0("LINC00211 vs SOFA (", label, ", n=", nrow(clin), ")"),
         x = "SOFA Score", y = "LINC00211 (log2CPM)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")
  plot_list[[label]] <- p
}

p_sofa <- ggarrange(plot_list[["Pediatric"]], plot_list[["Adult"]],
                    ncol = 2, common.legend = TRUE, legend = "right")
ggsave(file.path(out_dir, "Fig1_LINC00211_SOFA_scatter.pdf"),
       p_sofa, width = 10, height = 4.5)
ggsave(file.path(out_dir, "Fig1_LINC00211_SOFA_scatter.png"),
       p_sofa, width = 10, height = 4.5, dpi = 300)
cat("Saved: Fig1_LINC00211_SOFA_scatter\n")

# --- Fig 2: 共表达条形图 ---
coexpr_plot <- coexpr_all
coexpr_plot$gene <- factor(coexpr_plot$gene, levels = unique(coexpr_plot$gene))
coexpr_plot$category <- factor(coexpr_plot$category,
                               levels = c("Myeloid (positive)","PIRAT1 Target (mechanism)",
                                          "CellChat Ligand","T cell (negative control)",
                                          "IL16 (cross-module)"))

p_coexpr <- ggplot(coexpr_plot, aes(x = gene, y = spearman_rho, fill = category)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = c(0.3, -0.3), color = "grey60", linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~ cohort, ncol = 2) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Myeloid (positive)" = "#E64B35",
    "PIRAT1 Target (mechanism)" = "#3C5488",
    "CellChat Ligand" = "#F39B7F",
    "T cell (negative control)" = "#00A087",
    "IL16 (cross-module)" = "#7E6148"
  )) +
  labs(title = "LINC00211 Co-expression with Module Genes",
       x = NULL, y = "Spearman rho",
       fill = "Category") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "Fig2_LINC00211_coexpression_bar.pdf"),
       p_coexpr, width = 8, height = 7)
ggsave(file.path(out_dir, "Fig2_LINC00211_coexpression_bar.png"),
       p_coexpr, width = 8, height = 7, dpi = 300)
cat("Saved: Fig2_LINC00211_coexpression_bar\n")

# --- Fig 3: 偏相关 vs 原始相关 对比图 ---
compare_df <- merge(
  coexpr_all[, c("cohort","gene","category","spearman_rho","spearman_p")],
  partial_results[, c("cohort","gene","partial_rho","partial_p")],
  by = c("cohort","gene")
)
compare_df$category <- factor(compare_df$category,
    levels = c("Myeloid (positive)","PIRAT1 Target (mechanism)",
               "CellChat Ligand","T cell (negative control)","IL16 (cross-module)"))

compare_df$gene <- factor(compare_df$gene, levels = levels(coexpr_plot$gene))

compare_long <- compare_df %>%
  pivot_longer(cols = c(spearman_rho, partial_rho),
               names_to = "type", values_to = "rho") %>%
  mutate(type = recode(type, "spearman_rho" = "Raw", "partial_rho" = "SOFA-controlled"))

p_compare <- ggplot(compare_long, aes(x = gene, y = rho, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  facet_wrap(~ cohort, ncol = 2) +
  coord_flip() +
  scale_fill_manual(values = c("Raw" = "#7E6148", "SOFA-controlled" = "#B09C85")) +
  labs(title = "Raw vs SOFA-controlled Co-expression",
       x = NULL, y = "Spearman rho", fill = NULL) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "Fig3_raw_vs_partial.pdf"),
       p_compare, width = 8, height = 7)
ggsave(file.path(out_dir, "Fig3_raw_vs_partial.png"),
       p_compare, width = 8, height = 7, dpi = 300)
cat("Saved: Fig3_raw_vs_partial\n")

# --- Fig 4: 28天死亡 箱线图 ---
death_plot_list <- list()
for (label in c("Pediatric","Adult")) {
  expr <- if (label == "Pediatric") ped_expr else adult_expr
  clin <- if (label == "Pediatric") ped_clin else adult_clin
  lnc_vals <- as.numeric(expr["LINC00211", ])
  death <- clin$death_28

  valid <- !is.na(death)
  df <- data.frame(LINC00211 = lnc_vals[valid],
                   death = factor(death[valid], levels = c(0,1), labels = c("Alive","Deceased")))

  p <- ggboxplot(df, x = "death", y = "LINC00211", fill = "death",
                 add = "jitter", size = 0.5) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
    scale_fill_manual(values = c("Alive" = "#4DBBD5", "Deceased" = "#E64B35")) +
    labs(title = paste0("28-day Mortality (", label, ")"),
         x = NULL, y = "LINC00211 (log2CPM)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
  death_plot_list[[label]] <- p
}

p_death <- ggarrange(death_plot_list[["Pediatric"]], death_plot_list[["Adult"]],
                     ncol = 2)
ggsave(file.path(out_dir, "Fig4_LINC00211_death_boxplot.pdf"),
       p_death, width = 9, height = 4.5)
ggsave(file.path(out_dir, "Fig4_LINC00211_death_boxplot.png"),
       p_death, width = 9, height = 4.5, dpi = 300)
cat("Saved: Fig4_LINC00211_death_boxplot\n")
