rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

train <- "ba1_to_xbb"
test <- "xbb_to_end"

train <- "alpha_to_ba1"
test <- "ba1_to_end"

# train <- "ba5_to_xbb"
# test <- "pirola_to_end"

train <- "ba1_to_ba5"
test <- "ba5_to_xbb"

train <- "alpha_to_delta"
test <- "delta_to_ba1"

exp_list <- list(c("ba1_to_xbb", "xbb_to_end"),
                 c("alpha_to_ba1", "ba1_to_end"),
                 c("alpha_to_delta", "delta_to_ba1"),
                 c("ba1_to_ba5", "ba5_to_xbb"))

morsels <- foreach(exp = exp_list) %do% {
  # exp <- exp_list[[2]]
  train <- exp[1]
  test <- exp[2]
  res <- fread(str_glue("results/ML_out.020225/cross_dataset_results/train_{train}.test_{test}.results.csv"))
  res_linkage <- fread(str_glue("results/ML_out.020225/cross_dataset_results/train_{train}.test_{test}.partition.dprime_only.results.csv"))
  
  stat_df <- fread(str_glue("results/mutation_stats/{test}.dprime_stats.csv"))
  
  plot_df <- res %>%
    left_join(stat_df)
  
  plot_linkage_df <- res_linkage %>%
    left_join(stat_df)
  
  # rho <- signif(cor(plot_df$y_test, plot_df$y_pred, method = "spearman"), 2)
  r2 <- signif(1 - sum((plot_df$y_pred - plot_df$y_test)^2) / sum((plot_df$y_test - mean(plot_df$y_test))^2), 2)
  mae <- signif(mean(abs(plot_df$y_test - plot_df$y_pred)), 2)
  
  r2_linkage <- signif(1 - sum((plot_linkage_df$y_pred - plot_linkage_df$y_test)^2) / sum((plot_linkage_df$y_test - mean(plot_linkage_df$y_test))^2), 2)
  mae_linkage <- signif(mean(abs(plot_linkage_df$y_test - plot_linkage_df$y_pred)), 2)
  
  return(bind_rows(tibble(train = train, test = test, 
                          score = r2,
                          metric = "R2", 
                          linkage = F),
                   tibble(train = train, test = test, 
                          score = r2_linkage,
                          metric = "R2", 
                          linkage = T),
                   tibble(train = train, test = test, 
                          score = mae,
                          metric = "MAE", 
                          linkage = F),
                   tibble(train = train, test = test, 
                          score = mae_linkage,
                          metric = "MAE", 
                          linkage = T)))
}

bind_rows(morsels) %>%
  mutate(exp = str_glue("Train:{train}\nTest:{test}")) %>%
  ggplot(aes(x = exp, y = score, fill = linkage)) +
  geom_col(position = "dodge",
           color = "black") +
  facet_grid(rows = vars(metric)) +
  scale_fill_manual(values = c("grey", "pink")) +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        axis.title = element_text(face = "bold")) +
  geom_text(aes(label = score),
            position = position_dodge(width = 0.9),
            vjust = -0.1) +
  ylim(0, 0.8)

ggsave("results/prediction_out/linkage_versus_non_linkage_model_performance.cross_dataset.pdf", dpi = 600, height = 8, width = 12)
# plot_df %>%
#   arrange(desc(y_pred)) %>%
#   head(1000) %>% View()
#   summarise(sum(future_high)) %>% View()
# 
# count_df <- plot_df %>%
#   group_by(future_high) %>%
#   summarise(n = n())
# 
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
#   labs(x = "Future freq.>10%", y = "Predicted future fitness", 
#        title = str_glue("{dataset}"))
# 
# ggsave(str_glue("results/prediction_out/{dataset}.boxplot.pdf"), dpi = 600, width = 5, height = 4)
# 
# plot_df %>%
#   arrange(desc(y_pred)) %>%
#   head(50) %>%
#   summarise(n = sum(future_high))
# 
# linreg <- lm(y_pred ~ y_test,
#              data = plot_df)
# MASS::studres(linreg)
# 
# merged <- plot_df %>%
#   mutate(resid = MASS::studres(linreg))
# 
# merged %>%
#   group_by(future_high) %>%
#   summarise(median(resid))
# 
# wilcox.test(resid ~ future_high, data = merged)
# 
# merged %>%
#   ggplot(aes(x = future_high, y = resid)) +
#   geom_boxplot() +
#   geom_pwc()

monthly_df <- fread("results/allele_frequency_out/observed_out/after_ba1.monthly.csv")

mut_filt <- error_df %>%
  filter(future_high) 

plot_df %>%
  filter(!future_high) %>%
  ggplot(aes(x = y_test, y = y_pred)) +
  geom_bin2d(bins = 30) +
  paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", 
                                    direction = -1,
                                    trans = "log10") +
  geom_smooth(method = "lm",
              data = plot_df) +
  geom_point(aes(x = y_test, y = y_pred),
             fill = "turquoise3",
             color = "black",
             pch = 21,
             size = 3,
             data = plot_df %>% filter(future_high)) +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        axis.title = element_text(face = "bold"),
        legend.position = "top")

monthly_df %>%
  filter(mutation_name %in% mut_filt$mutation_name) %>%
  group_by(mutation_name) %>%
  summarise(n = sum(prop > 0.1)) %>%
  left_join(mut_filt) %>%
  ggplot(aes(x = n, y = diff)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm")

monthly_df %>%
  filter(mutation_name %in% mut_filt$mutation_name) %>%
  ggplot(aes(x = collection_month, y = prop, color = mutation_name)) +
  geom_line() +
  theme(legend.position = "none")
