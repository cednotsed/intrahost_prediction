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
file_list <- file_list[!grepl("raw|codon_freq", file_list)]
file_list

morsels <- foreach(file_name = file_list) %do% {
  d <- gsub(file_dir, "", file_name)
  d <- gsub(".results.csv", "", d)
  fread(file_name) %>%
    mutate(alias = d)
}

# Cross dataset
# file_dir2 <- "results/ML_out/cross_dataset_results/"
# file_list2 <- list.files(file_dir2, full.names = T)
# 
# morsels2 <- foreach(file_name = file_list2) %do% {
#   d <- gsub(file_dir2, "", file_name)
#   d <- gsub(".results.csv", "", d)
#   fread(file_name) %>%
#     mutate(alias = d)
# }
# cross_res <- bind_rows(morsels2) %>%
#   group_by(alias) %>%
#   summarise(R2 = cor(y_test, y_pred) ^ 2,
#             MAE = mean(abs(y_test - y_pred))) %>%
#   separate(alias, c(NA, "dataset", "experiment"), "\\.") %>%
#   pivot_longer(!c(dataset, experiment), names_to = "metric", values_to = "score")

merged <- bind_rows(morsels) %>%
  filter(!grepl("_to_|partition|spike|future|tonkin", alias)) %>%
  select(-fit_time, -score_time) %>%
  pivot_longer(!alias, names_to = "metric", values_to = "score") %>%
  mutate(metric = ifelse(metric == "test_r2", "R2", "MAE")) %>%
  separate(alias, c("dataset", "experiment1", "experiment2"), "\\.", remove = F) %>%
  mutate(dataset = factor(dataset, c("early", "alpha", "delta",
                                     "ba1", "ba5", "xbb",
                                     "pirola"))) %>%
  filter(is.na(experiment1)| experiment1 == "max_freq_only"| (experiment1 == "partition" & experiment2 == "dprime_only")) %>%
  mutate(experiment = case_when(is.na(experiment1) ~ "Full model",
                                experiment1 == "max_freq_only" ~ "Intrahost only")) %>%
  mutate(experiment = factor(experiment, c("Intrahost only", "Full model")))

merged %>%
  distinct(dataset, experiment1, experiment2)

merged %>%
  filter(metric == "R2") %>%
  ggplot(aes(x = dataset, y = score, fill = experiment)) +
  geom_boxplot(position = position_dodge(preserve = "single")) +
  # scale_fill_manual(values = c("turquoise", "grey", "#846D86FF", "#A7DBD8FF")) +
  # facet_grid(rows = vars(metric)) +
  theme_bw() +
  labs(x = "Dataset", y = "Nested cross-validation score (R2)", fill = "Model predictors") +
  theme(text = element_text(family = "sans"),
        axis.title = element_text(face = "bold")) +
geom_pwc()

ggsave("results/prediction_out/intrahost_versus_no_intrahost.pdf", dpi = 600, width = 6, height = 3)


merged %>%
  group_by(alias, metric, dataset) %>%
  summarise(mean_score = mean(score)) %>%
  filter(grepl("max", alias)|!grepl("partition", alias)) %>%
  filter(metric == "R2") %>%
  mutate(experiment = ifelse(grepl("max", alias), "no", "yes")) %>%
  select(dataset, experiment, mean_score) %>%
  pivot_wider(id_cols = dataset, names_from = experiment, values_from = mean_score) %>%
  mutate(fold_diff = (yes - no) / no * 100)



