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

exps <- list(c("alpha_to_ba1", "ba1_to_end"),
             c("alpha_to_delta", "delta_to_ba1"),
             c("ba1_to_ba5", "ba5_to_xbb"),
             c("ba1_to_xbb", "xbb_to_end"))

mut_meta <- fread("results/allele_frequency_out/observed_out/all_time.aggregate.csv") %>%
  distinct(mutation_name, protein_name, codon_number)

morsels <- foreach(exp = exps) %do% {
  train <- exp[1]
  test <- exp[2]
  
  res <- fread(str_glue("results/ML_out.020225/cross_dataset_results/train_{train}.test_{test}.results.csv"))
  stat_df <- fread(str_glue("results/mutation_stats/{test}.stats.csv"))
  
  res %>%
    left_join(stat_df) %>%
    mutate(train = train, 
           test = test)
}

plot_df <- bind_rows(morsels) %>%
  mutate(is_rbd = ifelse(protein_name == "Spike",
                         "Spike", "non-spike")) %>%
  mutate(abs_error = abs(y_pred - y_test)) %>%
  mutate(experiment = str_glue("train:{train}\ntest:{test}"))

plot_df %>%
  ggplot(aes(x = experiment, y = abs_error, fill = is_rbd)) + 
  geom_boxplot() +
  geom_pwc() 

findoutlier <- function(x) {
  return(x > 2 * median(x))
}

plot_df %>%
  group_by(experiment) %>%
  mutate(is_outlier = findoutlier(abs_error)) %>%
  ggplot(aes(x = experiment, y = abs_error, fill = is_rbd)) +
  geom_boxplot() +
  theme_bw() +
  geom_pwc()
  
ggsave("results/prediction_out/spike_vs_nonspike.results.pdf", dpi = 600, height = 5, width = 8)

test <- plot_df %>%
  group_by(experiment) %>%
  mutate(is_outlier = findoutlier(abs_error)) %>%
  filter(future_high) %>%
  group_by(experiment, is_outlier, is_rbd) %>%
  summarise(n = n()) %>%
  ungroup() 

exp_list <- deframe(test %>% distinct(experiment))

foreach(exp = exp_list) %do% {
  exp = exp_list[1]
  mat <- test %>%
    filter(experiment == exp) %>%
    select(-experiment) %>%
    pivot_wider(id_cols = is_outlier, names_from = is_rbd, values_from = n) %>%
    column_to_rownames("is_outlier")
  
  print(exp)
  print(fisher.test(mat))
}
