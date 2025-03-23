rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

meta <- fread("data/metadata/all_sra_metadata.csv")
fna <- readDNAStringSet("data/alignments/reassembled.gap_filtered.pango_filtered.aln")
file_dir <- "results/pipeline_out.020225/codon_out/"
file_list <- list.files(file_dir, full.names = T)

codon_df <- fread("results/allele_frequency_out/codon_out/missense_freq.filt.csv.gz")

test <- fread("results/allele_frequency_out/codon_out/delta.csv")

merged <- codon_df %>%
  left_join(meta %>% select(id = biosample, dataset))

# Distinct mutations per dataset
merged %>%
  group_by(dataset) %>%
  summarise(n = n_distinct(mutation_name))

# Median per library
mut_per_lib <- merged %>%
  group_by(id) %>%
  summarise(n = n_distinct(mutation_name)) %>%
  ungroup()

tibble(id = names(fna)) %>%
  left_join(mut_per_lib) %>%
  left_join(meta %>% select(id = biosample, dataset)) %>%
  mutate(n = replace_na(n, 0)) %>%
  arrange(n) %>%
  group_by(dataset) %>%
  summarise(median_n = median(n))

