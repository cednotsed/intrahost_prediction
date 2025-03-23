rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(paletteer)

fna <- readDNAStringSet("data/alignments/reassembled.gap_filtered.pango_filtered.aln")
meta <- fread("data/metadata/all_sra_metadata.csv") %>%
  filter(biosample %in% names(fna))

library_counts <- meta %>%
  group_by(alias) %>%
  summarise(n_biosamples = n_distinct(biosample))

file_dir <- "results/ML_out.020225/results_out/"
file_list <- list.files(file_dir, full.names = T)
file_list

morsels <- foreach(file_name = file_list) %do% {
  d <- gsub(file_dir, "", file_name)
  d <- gsub(".results.csv", "", d)
  fread(file_name) %>%
    mutate(alias = d)
}


merged <- bind_rows(morsels) %>%
  filter(!grepl("max|dprime|tonkin", alias))
merged %>% distinct(alias)

plot_df <- merged %>%
  select(alias, test_r2, test_mae) %>%
  separate(alias, c("dataset", "exp1"), "\\.") %>%
  filter(dataset %in% c("early", "alpha", "delta", "ba1", "ba5", "xbb", "pirola")) %>%
  filter(is.na(exp1)| exp1 == "spike") %>%
  mutate(exp_annot = ifelse(is.na(exp1), "All SAVs", "Spike only")) %>%
  select(-exp1) %>%
  pivot_longer(!c(dataset, exp_annot), names_to = "metric", values_to = "score")

plt1 <- plot_df %>%
  mutate(metric = ifelse(metric == "test_mae", "MAE", "r2")) %>%
  mutate(dataset = factor(dataset, c("early", "alpha", "delta", "ba1", "ba5", "xbb", "pirola"))) %>%
  ggplot(aes(x = dataset, y = score, fill = exp_annot)) +
  geom_boxplot() +
  facet_grid(rows = vars(metric), scale = "free") +
  labs(x = "Timeframe", y = "Model performance", 
       fill = "SAV subset") +
  geom_pwc() +
  theme_bw() 

## PREDICTION ERRORS
file_dir <- "results/ML_out.020225/within_dataset_results/"
file_list <- list.files(file_dir, full.names = T)

morsels2 <- foreach(file_name = file_list) %do% {
  d <- gsub(file_dir, "", file_name)
  d <- gsub(".results.csv", "", d)
  fread(file_name) %>%
    mutate(alias = d)
}

merged2 <- bind_rows(morsels2) %>% 
  filter(!grepl("max|dprime|tonkin", alias)) %>%
  filter(alias %in% c("early", "alpha", "delta", "ba1", "ba5", "xbb", "pirola")) %>%
  mutate(spike = grepl("Spike", mutation_name))

plt2 <- merged2 %>%
  mutate(alias = alias, factor(alias, c("early", "alpha", "delta", "ba1", "ba5", "xbb", "pirola"))) %>%
  mutate(spike = ifelse(spike, "Spike", "non-Spike")) %>%
  ggplot(aes(x = alias, y = abs(y_test - y_pred), fill = spike)) +
  geom_boxplot() +
  labs(x = "Timeframe", y = "Absolute prediction error", 
       fill = "SAV subset") +
  geom_pwc() +
  theme_classic()

plt1

ggsave("results/prediction_out/spike_vs_non_spike.pdf", dpi = 600, height = 5, width = 8)

