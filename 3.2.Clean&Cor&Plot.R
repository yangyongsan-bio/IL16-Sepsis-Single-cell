# 公共清洗函数：返回清洗后的数据框和剔除数量
clean_outliers <- function(df, cluster_name, clinical_var = "SOFA", 
                           max_outlier_pct = 0.3) {
  temp <- df %>% 
    select(y = all_of(cluster_name), x = all_of(clinical_var)) %>%
    drop_na()
  
  n_total <- nrow(temp)
  if (n_total < 5) return(list(clean_data = NULL, n_outliers = 0))
  
  fit <- lm(y ~ x, data = temp)
  cooksd <- cooks.distance(fit)
  stud_res <- rstudent(fit)
  
  cutoff_cook <- 4 / n_total
  cutoff_res <- 2
  outlier_idx <- which(cooksd > cutoff_cook | abs(stud_res) > cutoff_res)
  
  max_outliers <- floor(n_total * max_outlier_pct)
  if (length(outlier_idx) > max_outliers) {
    outlier_order <- order(cooksd, decreasing = TRUE)
    outlier_idx <- outlier_order[1:max_outliers]
  }
  
  temp_clean <- temp[-outlier_idx, , drop = FALSE]
  list(clean_data = temp_clean, n_outliers = length(outlier_idx))
}

# 函数1：自动清洗并计算 Spearman 相关（用于批量统计）
auto_clean_and_cor <- function(df, cluster_name, clinical_var = "SOFA", 
                               max_outlier_pct = 0.3) {
  result <- clean_outliers(df, cluster_name, clinical_var, max_outlier_pct)
  temp_clean <- result$clean_data
  n_total <- nrow(df %>% select(y = all_of(cluster_name), x = all_of(clinical_var)) %>% drop_na())
  n_outliers <- result$n_outliers
  n_clean <- ifelse(is.null(temp_clean), 0, nrow(temp_clean))
  
  # 原始（不清洗）的相关性
  temp_raw <- df %>% 
    select(y = all_of(cluster_name), x = all_of(clinical_var)) %>%
    drop_na()
  cor_raw <- tryCatch(
    cor.test(temp_raw$y, temp_raw$x, method = "spearman", exact = FALSE),
    error = function(e) NULL
  )
  
  # 清洗后的相关性
  cor_clean <- NULL
  if (!is.null(temp_clean) && nrow(temp_clean) >= 3) {
    cor_clean <- tryCatch(
      cor.test(temp_clean$y, temp_clean$x, method = "spearman", exact = FALSE),
      error = function(e) NULL
    )
  }
  
  tibble(
    subcluster = cluster_name,
    n_total = n_total,
    n_outliers_removed = n_outliers,
    n_clean = n_clean,
    rho_raw = ifelse(is.null(cor_raw), NA, cor_raw$estimate),
    p_raw = ifelse(is.null(cor_raw), NA, cor_raw$p.value),
    rho_clean = ifelse(is.null(cor_clean), NA, cor_clean$estimate),
    p_clean = ifelse(is.null(cor_clean), NA, cor_clean$p.value)
  )
}

# 函数2：获取清洗后的数据（用于绘图）
auto_clean_and_plot <- function(df, cluster_name, clinical_var = "SOFA", 
                             max_outlier_pct = 0.3) {
  result <- clean_outliers(df, cluster_name, clinical_var, max_outlier_pct)
  if (is.null(result$clean_data)) return(NULL)
  result$clean_data$subcluster <- cluster_name
  result$clean_data
}