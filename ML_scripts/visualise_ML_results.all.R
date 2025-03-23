rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(paletteer)

std <- function(x) sd(x)/sqrt(length(x))

fna <- readDNAStringSet("data/alignments/reassembled.gap_filtered.pango_filtered.aln")
meta <- fread("data/metadata/all_sra_metadata.csv") %>%
  filter(biosample %in% names(fna))

library_counts <- meta %>%
  group_by(alias) %>%
  summarise(n_biosamples = n_distinct(biosample))

file_dir <- "results/ML_out.020225/results_out/"
file_list <- list.files(file_dir, full.names = T)
file_list

within_morsels <- foreach(file_name = file_list) %do% {
  d <- gsub(file_dir, "", file_name)
  d <- gsub(".results.csv", "", d)
  fread(file_name) %>%
    mutate(alias = d)
}

within_df <- bind_rows(within_morsels) %>%
  filter(!grepl("max|partition|spike|tonkin", alias)) %>%
  select(alias, test_r2, test_mae) %>%
  group_by(alias) %>%
  summarise(se_mae = std(test_mae),
            test_mae = mean(test_mae),
            se_r2 = std(test_r2),
            test_r2 = mean(test_r2))

exp_list <- list(c("ba1_to_xbb", "xbb_to_end"),
                 c("alpha_to_ba1", "ba1_to_end"),
                 c("alpha_to_delta", "delta_to_ba1"),
                 c("ba1_to_ba5", "ba5_to_xbb"))

cross_morsels <- foreach(exp = exp_list) %do% {
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
  
  res <- bind_rows(train = train,
                   test = test,
                   mae = mae,
                   rho = rho,
                   r2 = r2) %>%
    mutate(exp = str_glue("Train: {train}\nTest: {test}"))
  
  return(res)
}

cross_df <- bind_rows(cross_morsels) %>%
  select(alias = exp, test_mae = mae, test_r2 = r2)

merged <- bind_rows(within_df %>%
                      mutate(type = "within"), 
                    cross_df %>%
                      mutate(type = "cross")) %>%
  pivot_longer(!c(alias, type), names_to = "metric", values_to = "score") %>%
  filter(!grepl("se_", metric)) %>%
  mutate(metric = ifelse(metric == "test_r2", "R2", "MAE")) %>%
  mutate(alias = factor(alias, c("early", "alpha", "delta",
                                 "ba1", "ba5", "xbb",
                                 "pirola", 
                                 "Train: alpha_to_ba1\nTest: ba1_to_end",
                                 "Train: alpha_to_delta\nTest: delta_to_ba1",
                                 "Train: ba1_to_ba5\nTest: ba5_to_xbb",
                                 "Train: ba1_to_xbb\nTest: xbb_to_end")))

se_df <- within_df %>%
  select(alias, se_mae, se_r2) %>%
  pivot_longer(!alias, names_to = "metric", values_to = "se") %>%
  mutate(metric = ifelse(metric == "se_r2", "R2", "MAE")) %>%
  mutate(type = "within")

plot_df <- merged %>% 
  left_join(se_df) %>%
  mutate(se = replace_na(se, 0)) %>%
  mutate(type = factor(type, c("within", "cross"))) %>%
  mutate(metric = factor(metric, c("R2", "MAE"))) %>%
  mutate(alias = factor(alias, c("early", "alpha", "delta",
                                 "ba1", "ba5", "xbb",
                                 "pirola", 
                                 "Train: alpha_to_ba1\nTest: ba1_to_end",
                                 "Train: alpha_to_delta\nTest: delta_to_ba1",
                                 "Train: ba1_to_ba5\nTest: ba5_to_xbb",
                                 "Train: ba1_to_xbb\nTest: xbb_to_end")))

plot_df %>%
  ggplot(aes(x = alias, y = score, fill = metric)) +
  geom_col(color = "black") +
  geom_errorbar(aes(x = alias, 
                    ymin = score-se, 
                    ymax = score+se), width=0.2, colour="black") +
  facet_grid(rows = vars(metric),
             cols = vars(type), 
             space = "free",
             scales = "free") +
  theme_bw() 

ggsave(str_glue("results/ML_out.020225/all_results.pdf"), dpi = 600, width = 10, height = 7)

merged %>%
  filter(type == "within") %>% View()
  group_by(metric) %>%
  summarise(min = min(score),
            max = max(score))
# merged %>%
#   filter(experiment %in% c("Spike", "None", "Dprime only", "Partition + dprime", "Partition only")) %>%
#   filter(metric == "R2") %>%
#   ggplot(aes(x = dataset, y = score, fill = experiment)) +
#   geom_boxplot(position = position_dodge(preserve = "single")) +
#   # scale_fill_manual(values = c("turquoise", "grey", "#846D86FF", "#A7DBD8FF")) +
#   # facet_grid(rows = vars(metric)) +
#   theme_bw() +
#   labs(x = "Dataset", y = "Nested cross-validation score") +
#   theme(text = element_text(family = "sans"),
#         axis.title = element_text(face = "bold"))
#   # geom_pwc()
# 
# #   ggplot(aes(x = log10(global_n + 1))) +
# #   geom_histogram()
# 
# merged %>%
#   group_by(dataset, experiment, metric) %>%
#   summarise(mean_score = mean(score)) %>%
#   filter(experiment == "Spike") %>%
#   filter(metric == "MAE") %>% View()
# 
# merged %>%
#   filter(experiment == "Spike") %>% View()