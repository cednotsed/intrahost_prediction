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

merged <- bind_rows(morsels) %>%
  filter(!grepl("_to_", alias)) %>%
  select(-fit_time, -score_time) %>%
  pivot_longer(!alias, names_to = "metric", values_to = "score") %>%
  mutate(metric = ifelse(metric == "test_r2", "R2", "MAE")) %>%
  separate(alias, c("dataset"), "\\.", remove = F) %>% 
  rowwise() %>%
  mutate(experiment = gsub(str_glue("{dataset}"), "", alias)) %>%
  mutate(experiment = gsub("\\.", "", experiment))

merged %>%
  distinct(experiment)
  # # distinct(experiment1, experiment2)
  mutate(experiment = case_when(experiment1 == "dprime" ~ "Dprime + corr.",
                                experiment1 == "dprime_only" ~ "Dprime only",
                                experiment1 == "partition" & experiment2 == "dprime_only" ~ "Partition + dprime",
                                experiment1 == "partition" & is.na(experiment2) ~ "Partition only",
                                experiment1 == "linkage" ~ "Corr. only",
                                experiment1 == "spike" ~ "Spike",
                                experiment1 == "prior" & is.na(experiment2) ~ "Prior",
                                is.na(experiment1) & is.na(experiment2) ~ "None",
                                experiment1 == "prior" & experiment2 == "dprime_only" ~ "Prior + Dprime")) %>%
  mutate(dataset = factor(dataset, c("early", "alpha", "delta",
                                     "ba1", "ba5", "xbb",
                                     "pirola")))

merged %>%
  distinct(experiment1, experiment2, experiment)

merged %>%
  filter(experiment %in% c("Spike", "None", "Dprime only", "Partition + dprime", "Partition only")) %>%
  filter(metric == "R2") %>%
  ggplot(aes(x = dataset, y = score, fill = experiment)) +
  geom_boxplot(position = position_dodge(preserve = "single")) +
  # scale_fill_manual(values = c("turquoise", "grey", "#846D86FF", "#A7DBD8FF")) +
  # facet_grid(rows = vars(metric)) +
  theme_bw() +
  labs(x = "Dataset", y = "Nested cross-validation score") +
  theme(text = element_text(family = "sans"),
        axis.title = element_text(face = "bold"))
  # geom_pwc()

#   ggplot(aes(x = log10(global_n + 1))) +
#   geom_histogram()

merged %>%
  group_by(dataset, experiment, metric) %>%
  summarise(mean_score = mean(score)) %>%
  filter(experiment == "Spike") %>%
  filter(metric == "MAE") %>% View()

merged %>%
  filter(experiment == "Spike") %>% View()