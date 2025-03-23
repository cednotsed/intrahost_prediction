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

merged <- bind_rows(morsels) %>%
  mutate(is_rbd = ifelse(protein_name == "Spike",
                         "Spike", "non-spike")) %>%
  mutate(error = y_pred - y_test,
         abs_error = abs(y_pred - y_test)) %>%
  mutate(experiment = str_glue("train:{train}\ntest:{test}"))

morsels <- foreach(exp = unique(merged$experiment)) %do% {
  temp <- merged %>%
    filter(experiment == exp) 
  
  linreg <- lm(y_pred ~ y_test,
               data = temp)
  
  temp_parsed <- temp %>%
    mutate(resid = MASS::studres(linreg))
  
  return(temp_parsed)
}

bind_rows(morsels) %>%
  ggplot(aes(x = experiment, y = resid, fill = future_high)) +
  geom_boxplot() +
  geom_pwc() +
  theme_bw() +
  scale_fill_manual(values = c("indianred", "royalblue")) +
  labs(x = "Cross-dataset model", y = "Studentised residuals") +
  theme(legend.position = "top")

ggsave("results/prediction_out/fit_vs_nonfit.studres.pdf", dpi = 600, width = 8, height = 5)

