rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

meta <- fread("data/metadata/all_sra_metadata.csv")

timeframes <- meta %>%
  group_by(dataset, alias) %>%
  summarise(start = min(collection_month),
            end = max(collection_month)) %>%
  mutate(start = start %m-% months(1),
         end = end %m+% months(1)) %>%
  mutate(end_month = as.numeric(format(end, "%m")),
         end_year = as.numeric(format(end, "%Y")))

monthly_df <- fread("results/allele_frequency_out/observed_out/all_monthly_frequencies.csv") %>%
  mutate(collection_date = as.Date(str_glue("{collection_month}-01")))

exp_list <- list(c("ba1_to_xbb", "xbb_to_end"),
                 c("alpha_to_ba1", "ba1_to_end"),
                 c("alpha_to_delta", "delta_to_ba1"),
                 c("ba1_to_ba5", "ba5_to_xbb"),
                 c("tonkin_to_delta", "ba1_to_end.US_only"),
                 c("early_to_delta.non_US", "ba1_to_end.US_only"))

morsels <- foreach(exp = exp_list) %do% {
  train <- exp[1]
  test <- exp[2]
  
  timeframe <- timeframes %>%
    filter(alias == str_split(test, "\\_")[[1]][1])
  
  monthly_filt <- monthly_df %>%
    filter(collection_month <= timeframe$end)
  
  res <- fread(str_glue("results/ML_out.020225/cross_dataset_results/train_{train}.test_{test}.results.csv"))
  
  stat_df <- fread(str_glue("results/mutation_stats/{test}.stats.csv"))
  
  res %>%
    left_join(stat_df) %>%
    mutate(type = ifelse(mutation_name %in% monthly_filt$mutation_name, "Observed", "Unobserved")) %>%
    mutate(exp = str_glue("Train: {train}\nTest: {test}"))
  
}

merged <- bind_rows(morsels) %>%
  mutate(prediction_error = abs(y_test - y_pred))

merged %>%
  ggplot(aes(x = exp, y = prediction_error, fill = type)) +
  geom_boxplot() +
  theme_bw() +
  geom_pwc() +
  labs(x = "Cross-dataset model", y = "Abs. prediction error", fill = "SAV type")

ggsave("results/prediction_out/observed_vs_unobserved.pdf", dpi = 600, width = 12, height = 5)
