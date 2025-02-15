rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)

fna <- readDNAStringSet("data/alignments/reassembled.gap_filtered.pango_filtered.aln")
df <- fread("data/metadata/all_sra_metadata.csv") %>%
  select(biosample, bases, dataset) %>%
  filter(biosample %in% names(fna)) %>%
  filter(dataset == "Early")

tonkin_dedup <- fread("data/external_datasets/tonkin/tonkin.dedup.txt")

tonkin_df <- fread("data/external_datasets/tonkin/tonkin.metadata.csv") %>%
  filter(sample_name %in% tonkin_dedup$sample_name) %>%
  mutate(dataset = "Tonkin") %>%
  select(sample_name, bases, dataset)

bind_rows(tonkin_df, df) %>%
  ggplot(aes(x = dataset, y = log10(bases))) +
  geom_boxplot()
  group_by(dataset) %>%
  summarise(median = median(bases))
