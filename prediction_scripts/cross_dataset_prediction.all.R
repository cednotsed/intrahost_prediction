rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

train <- "alpha_to_ba1"
test <- "ba1_to_end"

exp_list <- list(c("ba1_to_xbb", "xbb_to_end"),
                 c("alpha_to_ba1", "ba1_to_end"),
                 c("alpha_to_delta", "delta_to_ba1"),
                 c("ba1_to_ba5", "ba5_to_xbb"))

morsels <- foreach(exp = exp_list) %do% {
  # exp <- exp_list[[2]]
  train <- exp[1]
  test <- exp[2]
  
  res <- fread(str_glue("results/ML_out.020225/cross_dataset_results/train_{train}.test_{test}.results.csv"))
  stat_df <- fread(str_glue("results/mutation_stats/{test}.stats.csv"))
  
  plot_df <- res %>%
    left_join(stat_df)
  
  mae <- signif(mean(abs(plot_df$y_test - plot_df$y_pred)), 2)
  rho <- signif(cor(plot_df$y_test, plot_df$y_pred, method = "spearman"), 2)
  r2 <- signif(1 - sum((plot_df$y_pred - plot_df$y_test)^2) / sum((plot_df$y_test - mean(plot_df$y_test))^2), 2)
  
  # plot_df %>%
  #   filter(!future_high) %>%
  #   ggplot(aes(x = y_test, y = y_pred)) +
  #   geom_bin2d(bins = 30) +
  #   paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", 
  #                                     direction = -1,
  #                                     trans = "log10") +
  #   geom_smooth(method = "lm",
  #               data = plot_df) +
  #   geom_point(aes(x = y_test, y = y_pred),
  #              fill = "turquoise3",
  #              color = "black",
  #              pch = 21,
  #              size = 3,
  #              data = plot_df %>% filter(future_high)) +
  #   theme_bw() +
  #   theme(text = element_text(family = "sans"),
  #         axis.title = element_text(face = "bold"),
  #         # legend.position = "top") +
  #         legend.position = "none") +
  #   ylim(-0.6, 6) + 
  #   labs(x = "Observed fitness", y = "Predicted fitness",
  #        title = str_glue("{test}: rho={rho}, r2={r2}, MAE={mae}")) 
  
  # ggsave(str_glue("results/prediction_out/{train}.{test}.cross_data_results.pdf"), dpi = 600, width = 5, height = 3.5)
  
  # plot_df %>%
  #   arrange(desc(y_pred)) %>%
  #   head(1000) %>% View()
  #   summarise(sum(future_high)) %>% View()
  
  count_df <- plot_df %>%
    group_by(future_high) %>%
    summarise(n = n())
  
  # plot_df %>%
  #   arrange(desc(y_pred)) %>% 
  #   ggplot(aes(x = future_high, y = y_pred, fill = future_high)) +
  #   geom_boxplot() +
  #   scale_fill_manual(values = c("indianred", "turquoise3")) +
  #   geom_text(aes(x = future_high, y = max(plot_df$y_pred) + 0.5, label = str_glue("n={n}")),
  #             data = count_df) +
  #   geom_pwc() +
  #   theme_bw() +
  #   theme(text = element_text(family = "sans"),
  #         legend.position = "none",
  #         axis.title = element_text(face = "bold")) +
  #   labs(x = "Future freq.>10%", y = "Predicted future fitness")
  
  res <- bind_rows(train = train,
            test = test,
            mae = mae,
            rho = rho,
            r2 = r2) %>%
    mutate(exp = str_glue("Train: {train}\nTest: {test}"))
  
  return(res)
}

bind_rows(morsels) %>%
  ggplot(aes(x = exp, y = r2)) +

plot_df %>%
  ggplot(aes(x = alias))
ggsave(str_glue("results/prediction_out/{dataset}.boxplot.pdf"), dpi = 600, width = 5, height = 4)

plot_df %>%
  arrange(desc(y_pred)) %>%
  head(50) %>%
  summarise(n = sum(future_high))

linreg <- lm(y_pred ~ y_test,
             data = plot_df)
MASS::studres(linreg)

merged <- plot_df %>%
  mutate(resid = MASS::studres(linreg))

merged %>%
  group_by(future_high) %>%
  summarise(median(resid))

wilcox.test(resid ~ future_high, data = merged)

merged %>%
  ggplot(aes(x = future_high, y = resid)) +
  geom_boxplot() +
  geom_pwc()
