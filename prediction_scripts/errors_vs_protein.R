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
  mutate(abs_error = abs(y_pred - y_test)) %>%
  mutate(experiment = str_glue("train:{train}\ntest:{test}")) %>%
  mutate(region = case_when(grepl("NSP", protein_name) ~ "ORF1ab",
                            grepl("NS7", protein_name) ~ "NS7",
                            TRUE ~ protein_name)) %>%
  filter(protein_name != "NS10", protein_name != "NSP11") %>%
  mutate(region = factor(region, c("ORF1ab", "Spike", "NS3",
                                   "E", "M", "NS6", 
                                   "NS7", "NS8",
                                   "N"))) 

plot_df %>%
  ggplot(aes(x = region, y = abs_error, fill = region)) + 
  geom_boxplot(outliers = F) +
  facet_grid(rows = vars(experiment)) +
  theme_bw() +
  labs(x = "Region", y = "Abs. prediction error") +
  theme(legend.position = "none")

ggsave("results/prediction_out/errors_by_protein.pdf", dpi = 600, height = 6, width = 8)

