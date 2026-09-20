
suppressMessages(library(data.table))

# ---------------- 配色 ---------------- 
col_rep  <- "#C0392B"   # REPLICATED
col_cons <- "#E67E22"   # direction consistent (n.s.)
col_exp  <- "#7F8C8D"   # exploratory
col_not  <- "#95A5A6"   # not consistent

vorder <- c("REPLICATED" = 1, "direction consistent (n.s.)" = 2,
            "exploratory" = 3, "not consistent" = 4,
            "not assessable (no probe on platform)" = 5)

verdict_cols <- function(v) {
  ifelse(v == "REPLICATED", col_rep,
  ifelse(v == "direction consistent (n.s.)", col_cons,
  ifelse(v == "exploratory", col_exp, col_not)))
}

# ---------------- 条形图面板 ----------------  (genes 第一个在最上方)
draw_bar <- function(genes, delta, rb, pv, cols, xlim, xticks, main, xlab, note, leg_lab, leg_col) {
  n  <- length(genes)
  y  <- seq(n, 1)                                   
  bh <- 0.34                                        
  cexy <- if (n > 16) 0.62 else if (n > 10) 0.72 else 0.78
  par(mar = c(6.2, max(6.8, max(nchar(genes)) * 0.52 + 1.2), 4.6, 0.8))
  plot.new()
  plot.window(xlim = xlim, ylim = c(0.4, n + 2.3))  
  abline(v = 0, col = "#9AA5B1", lty = 2, lwd = 1.1)
  off <- 0.022 * diff(xlim)
  for (i in seq_len(n)) {
    rect(min(0, delta[i]), y[i] - bh, max(0, delta[i]), y[i] + bh,
         col = cols[i], border = NA)
    lab <- if (is.na(pv[i])) sprintf("r = %+.3f", rb[i])
           else sprintf("r = %+.3f %s", rb[i], if (pv[i] < 0.05) "*" else "ns")
    if (delta[i] >= 0) text(delta[i] + off, y[i], lab, adj = 0, cex = 0.66, col = "#2C3E50")
    else               text(delta[i] - off, y[i], lab, adj = 1, cex = 0.66, col = "#2C3E50")
  }
  axis(1, at = xticks, cex.axis = 0.78, mgp = c(1.8, 0.4, 0))
  axis(2, at = y, labels = genes, las = 1, cex.axis = cexy)
  mtext(xlab, side = 1, line = 3.0, cex = 0.78)
  mtext(note, side = 1, line = 4.3, cex = 0.64, adj = 0, col = "#5F5E5A")
  title(main = main, adj = 0, cex.main = 0.82, font.main = 1, line = 1.0)
  legend("top", legend = leg_lab, pch = 15, col = leg_col,
         cex = 0.68, bty = "n", ncol = 2, x.intersp = 0.8, y.intersp = 0.9)
}

# ---------------- 评分箱线图面板 ----------------
draw_score_box <- function(sc, death, xnames, main, ylab, seed = 42) {
  surv <- sc[death == 0]; nonsurv <- sc[death == 1]
  par(mar = c(4.4, 4.0, 3.4, 0.6))
  ylim <- range(c(sc, 0)) + c(-0.5, 0.5)
  boxplot(list(nonsurv, surv), names = xnames,
          col = c("#F5B7B1", "#AED6F1"), outline = FALSE, boxwex = 0.5,
          ylim = ylim, axes = FALSE, border = "#34495E",
          medcol = "#333333", medlwd = 1.6, staplewex = 0)
  abline(h = 0, col = "#9AA5B1", lty = 2, lwd = 1)
  axis(2, cex.axis = 0.78, mgp = c(1.8, 0.4, 0), las = 1)
  axis(1, at = 1:2, labels = xnames, cex.axis = 0.78, mgp = c(1.8, 0.55, 0))
  set.seed(seed)
  for (i in 1:2) {
    g  <- list(nonsurv, surv)[[i]]
    xs <- i + runif(length(g), -0.13, 0.13)
    points(xs, g, pch = 16, cex = 0.5,
           col = adjustcolor(c("#C0392B", "#2E86C1")[i], alpha.f = 0.4))
  }
  mtext(ylab, side = 2, line = 2.6, cex = 0.8)
  title(main = main, adj = 0, cex.main = 0.82, font.main = 1, line = 1.0)
}

# 一、成人队列 GSE65682  (n=479; 114 non-survivors / 365 survivors)
resA <- fread("GSE65682_validation/GSE65682_validation_results_R.csv")
# 仅纳入可评估的正式基因 (ND4L / LINC00211 平台无探针 -> 不进面板)
formalA <- resA[!is.na(primary_dir) & !is.na(delta_z)]
formalA <- formalA[order(vorder[verdict], p_wilcoxon)]
scoreA  <- fread("GSE65682_validation/GSE65682_score_data_R.csv")
deathA  <- scoreA$death; scA <- scoreA$score
nA      <- nrow(formalA)

pA <- wilcox.test(scA[deathA == 1], scA[deathA == 0],
                  correct = FALSE)$p.value        # two-sided

plot_adult <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.6, 1))
  draw_bar(formalA$gene, formalA$delta_z, formalA$rank_biserial_r, formalA$p_wilcoxon,
           verdict_cols(formalA$verdict),
           xlim   = c(-0.56, 1.32), xticks = seq(-0.4, 0.8, by = 0.2),
           main   = "A  Module genes and 28-day mortality\nGSE65682 cohort (n=479; 114 non-survivors / 365 survivors)",
           xlab   = "Standardized mean difference (non-survivors - survivors), z units",
           note   = "* p < 0.05; ns = not significant;  r = rank-biserial r (negative: higher in non-survivors);  ND4L / LINC00211 not assessable",
           leg_lab = c("Replicated (p<0.05)", "Direction consistent (n.s.)"),
           leg_col = c(col_rep, col_cons))
  draw_score_box(scA, deathA, c("Non-survivors (n=114)", "Survivors (n=365)"),
                 main = sprintf("B  Direction-weighted polygenic score vs mortality\nMann-Whitney p = %s (two-sided)",
                                sprintf("%.2e", pA)),
                 ylab = sprintf("Direction-weighted polygenic score (weighted z, %d genes)", nA))
}

# 二、儿童队列 GSE66099  (n=108; 13 non-survivors / 95 survivors)
resP <- fread("GSE66099_validation/GSE66099_validation_results_R.csv")
resP <- resP[!is.na(delta_z)]                       
resP <- resP[order(vorder[verdict], p_wilcoxon)]    
scoreP <- fread("GSE66099_validation/GSE66099_score_data_R.csv")
deathP <- scoreP$death; scP <- scoreP$score
nP     <- nrow(resP)                                
nP_score <- sum(!is.na(resP$primary_dir) & resP$primary_dir != 0)   

pP <- wilcox.test(scP[deathP == 1], scP[deathP == 0], correct = FALSE)$p.value

plot_pediatric <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.6, 1))
  draw_bar(resP$gene, resP$delta_z, resP$rank_biserial_r, resP$p_wilcoxon,
           verdict_cols(resP$verdict),
           xlim   = c(-0.56, 1.32), xticks = seq(-0.4, 0.8, by = 0.2),
           main   = "A  Module genes and 28-day mortality\nGSE66099 cohort (n=108; 13 non-survivors / 95 survivors)",
           xlab   = "Standardized mean difference (non-survivors - survivors), z units",
           note   = "* p < 0.05; ns = not significant;  r = rank-biserial r (negative: higher in non-survivors);  LINC00211 not assessable",
           leg_lab = c("Replicated (p<0.05)", "Direction consistent (n.s.)", "Not consistent"),
           leg_col = c(col_rep, col_cons, col_not))
  draw_score_box(scP, deathP, c("Non-survivors (n=13)", "Survivors (n=95)"),
                 main = sprintf("B  Direction-weighted polygenic score vs mortality\nMann-Whitney p = %s (two-sided)",
                                sprintf("%.2e", pP)),
                 ylab = sprintf("Direction-weighted polygenic score (weighted z, %d genes)", nP_score))
}


# 三、输出 PNG + PDF 
out_specs <- list(
  list(tag = "GSE65682_validation/GSE65682_validation_figure_R",
       plot = plot_adult,      h = max(5.2, 1.9 + 0.40 * nA)),
  list(tag = "GSE66099_validation/GSE66099_validation_figure_R",
       plot = plot_pediatric,  h = max(5.2, 1.9 + 0.40 * nP))
)
for (sp in out_specs) {
  png(paste0(sp$tag, ".png"), width = 11, height = sp$h, units = "in", res = 300)
  sp$plot(); dev.off()
  pdf(paste0(sp$tag, ".pdf"), width = 11, height = sp$h)
  sp$plot(); dev.off()
}
cat(sprintf("figures saved (adult %d genes, %.1f in; pediatric %d genes, %.1f in)\n",
            nA, out_specs[[1]]$h, nP, out_specs[[2]]$h))
