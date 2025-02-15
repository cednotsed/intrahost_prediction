rm(list = ls())
setwd("C:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

fna <- readDNAStringSet("data/alignments/reassembled.gap_filtered.pango_filtered.aln")
file_dir <- "results/pipeline_out.020225/stats_out/"
file_list <- list.files(file_dir, "coverage", full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  fread(file_name) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(biosample = file_name)
}

merged <- bind_rows(morsels) %>%
  mutate(biosample = gsub(".coverage.txt", "", biosample)) %>%
  mutate(biosample = gsub(file_dir, "", biosample)) %>%
  mutate(meanbaseq = as.numeric(meanbaseq)) %>%
  mutate(meanmapq = as.numeric(meanmapq)) %>%
  filter(biosample %in% names(fna))

nrow(merged) == length(fna)

merged %>%
  ggplot(aes(x = as.numeric(meanbaseq))) +
  geom_histogram()

merged %>%
  ggplot(aes(x = as.numeric(meanmapq))) +
  geom_histogram()

min_score <- deframe(merged %>%
  summarise(min(meanbaseq)))

10 ^ (-min_score / 10)
