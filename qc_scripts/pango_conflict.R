rm(list = ls())
setwd("C:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

fna <- readDNAStringSet("data/alignments/reassembled.gap_filtered.pango_filtered.aln")
meta <- fread("data/metadata/all_sra_metadata.csv")

pango_df <- fread("data/metadata/reassembled.gap_filtered.pango_lineage.csv") %>%
  dplyr::rename(biosample = taxon) %>%
  filter(biosample %in% names(fna)) %>%
  left_join(meta %>% select(biosample, dataset))

pango_df %>%
  filter(scorpio_conflict > 0) %>% View()
  ggplot(aes(x = scorpio_conflict)) +
  geom_histogram()
  filter(scorpio_conflict > 0)
