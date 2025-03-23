rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

# train <- "ba1_to_xbb"
# test <- "xbb_to_end"
# 
# train <- "alpha_to_ba1"
# test <- "ba1_to_end"
# 
# # train <- "ba5_to_xbb"
# # test <- "pirola_to_end"
# 
# train <- "ba1_to_ba5"
# test <- "ba5_to_xbb"
# 
# train <- "alpha_to_delta"
# test <- "delta_to_ba1"

exp_list <- list(c("ba1_to_xbb", "xbb_to_end"),
                 c("alpha_to_ba1", "ba1_to_end"),
                 c("alpha_to_delta", "delta_to_ba1"),
                 c("ba1_to_ba5", "ba5_to_xbb"),
                 c("tonkin_to_delta", "ba1_to_end.US_only"),
                 c("early_to_delta.non_US", "ba1_to_end.US_only"))

morsels <- foreach(exp = exp_list) %do% {
  # exp <- exp_list[[1]]
  train <- exp[1]
  test <- exp[2]
  res <- fread(str_glue("results/ML_out.020225/cross_dataset_results/train_{train}.test_{test}.results.csv"))
  
  stat_df <- fread(str_glue("results/mutation_stats/{test}.stats.csv"))
  
  plot_df <- res %>%
    left_join(stat_df)
  
  rho <- signif(cor(plot_df$y_test, plot_df$y_pred, method = "spearman"), 2)
  r2 <- signif(1 - sum((plot_df$y_pred - plot_df$y_test)^2) / sum((plot_df$y_test - mean(plot_df$y_test))^2), 2)
  mae <- signif(mean(abs(plot_df$y_test - plot_df$y_pred)), 2)
  
  tibble(train = train, test = test,
         rho = rho,
         r2 = r2,
         mae = mae)
}

bind_rows(morsels)
