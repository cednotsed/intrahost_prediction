rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

d <- "xbb"
datasets <- c("alpha", "delta", "ba1", "ba5", "xbb", "pirola")
# datasets <- c("alpha", "delta", "ba1", "xbb", "pirola")

morsels <- foreach(d = datasets) %do% {
  res <- fread(str_glue("results/ML_out.020225/results_out/{d}.results.csv"))
  res_linkage <- fread(str_glue("results/ML_out.020225/results_out/{d}.partition.dprime_only.results.csv"))
  
  stat_df <- fread(str_glue("results/mutation_stats/{d}.dprime_stats.csv"))
  
  bind_rows(res %>% mutate(dataset = d, type = "No linkage"),
            res_linkage %>% mutate(dataset = d, type = "Linkage"))
}

plot_df <- bind_rows(morsels) %>%
  mutate(dataset = factor(dataset, c("alpha", "delta", "ba1", "ba5", "xbb", "pirola"))) %>%
  select(dataset, type, test_r2, test_mae) %>%
  pivot_longer(!c(dataset, type), names_to = "metric", values_to = "score") %>%
  mutate(metric = ifelse(metric == "test_r2", "R2", "MAE")) %>%
  mutate(metric = factor(metric, c("R2", "MAE")))

plot_df %>%
  filter(metric == "R2") %>% 
  ggplot(aes(x = dataset, y = score, fill = type)) +
  geom_boxplot() +
  geom_pwc() +
  scale_fill_manual(values = c("pink", "grey")) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(family = "sans"),
        axis.title = element_text(face = "bold")) +
  labs(x = "Dataset")

ggsave(str_glue("results/prediction_out/linkage_vs_no_linkage.pdf"), dpi = 600, height = 3, width = 4.5)

plot_df %>%
  filter(metric == "R2") %>%
  filter(type == "Linkage") %>%
  group_by(dataset, type) %>%
  summarise(mean = mean(score)) %>% View()
# count_df <- plot_df %>%
#   group_by(alias, future_high) %>%
#   summarise(n = n())
#   
# plot_df %>%
#   # mutate(diff = link_error - error) %>%
#   ggplot(aes(x = alias, y = diff, fill = future_high)) +
#   geom_boxplot() +
#   geom_text(aes(y = 2, label = str_glue("{n}")),
#             data = count_df,
#             position = position_dodge(width = .9),
#             size = 2) +
#   geom_pwc(label.size = 2) +
#   theme_bw() +
#   scale_fill_manual(values = c("indianred", "royalblue")) +
#   theme(text = element_text(family = "sans"),
#         axis.title = element_text(face = "bold")) +
#   labs(x = "Dataset", y = "Change in prediction error")

# fwrite(plot_df, "results/prediction_out/prediction_errors.csv")
# error_df %>%
#   filter(max_dprime != -1) %>%
#   ggplot(aes(x = log10(global_n + 1), y = diff, color = max_dprime)) +
#   geom_point() +
#   geom_smooth(color = "#E5AD4FFF") +
#   theme_bw() +
#   theme(text = element_text(family = "sans"),
#         # legend.position = "none",
#         axis.title = element_text(face = "bold")) +
#   labs(x = "Log10(future fitness + 1)", y = "Change in prediction error", title = "")
# ggplot(aes(x = future_high, y = abs(y_test - y_pred))) +
# scale_fill_manual(values = c("indianred", "turquoise3")) +
# geom_boxplot(aes(fill = linkage)) +
# geom_text(aes(label = str_glue("n={n}"), y = -0.5),
#           data = count_df) +
# geom_pwc(aes(fill = linkage)) +
# theme_bw() +
#   geom_text(aes(label = str_glue("n={n}"), y = -0.5),
#             data = count_df) +
# theme(text = element_text(family = "sans"),
#       # legend.position = "none",
#       axis.title = element_text(face = "bold")) +
# labs(x = "Future freq.>10%", y = "Absolute prediction error", fill = "Linkage")
# 
# ggsave(str_glue("results/prediction_out/linkage_vs_no_linkage/prediction_errors.{d}.pdf"), dpi = 600, height = 3, width = 3.5)
# 
# return(tibble(train = train, test = test, 
#               r2 = r2, r2_linkage = r2_linkage, 
#               mae = mae, mae_linkage = mae_linkage))

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
# 
# monthly_df <- fread("results/allele_frequency_out/observed_out/after_ba1.monthly.csv")
# 
# mut_filt <- error_df %>%
#   filter(future_high) 
# 
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
#         legend.position = "top")
# 
# monthly_df %>%
#   filter(mutation_name %in% mut_filt$mutation_name) %>%
#   group_by(mutation_name) %>%
#   summarise(n = sum(prop > 0.1)) %>%
#   left_join(mut_filt) %>%
#   ggplot(aes(x = n, y = diff)) +
#   geom_point() +
#   theme_bw() +
#   geom_smooth(method = "lm")
# 
# monthly_df %>%
#   filter(mutation_name %in% mut_filt$mutation_name) %>%
#   ggplot(aes(x = collection_month, y = prop, color = mutation_name)) +
#   geom_line() +
#   theme(legend.position = "none")
