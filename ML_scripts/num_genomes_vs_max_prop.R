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
  mutate(experiment = gsub("\\.", "", experiment)) %>%
  mutate(experiment = case_when(experiment == "" ~ "No. genomes",
                                experiment == "max_future_prop" ~ "Max. monthly prop."))


merged %>%
  mutate(dataset = factor(dataset, c("early", "alpha", "delta", "ba1", "ba5", "xbb", "pirola"))) %>%
  filter(dataset != "tonkin") %>%
  filter(!is.na(experiment)) %>%
  filter(metric == "R2") %>%
  ggplot(aes(x = dataset, y = score, fill = experiment)) +
  geom_boxplot(position = position_dodge(preserve = "single")) +
  theme_bw() +
  labs(x = "Dataset", y = "Nested cross-validation score", fill = "Future fitness metric") +
  theme(text = element_text(family = "sans"),
        axis.title = element_text(face = "bold")) +
  geom_pwc()

ggsave("results/prediction_out/num_genomes_vs_max_prop.pdf", dpi = 600, height = 5, width = 8)


