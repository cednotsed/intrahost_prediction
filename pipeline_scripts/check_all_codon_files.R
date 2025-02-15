rm(list = ls())
setwd("C:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

meta <- fread("data/metadata/all_sra_metadata.csv")
fna <- readDNAStringSet("data/alignments/reassembled.gap_filtered.pango_filtered.aln")
file_dir <- "results/pipeline_out.020225/codon_out/"
file_list <- list.files(file_dir, full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  temp <- fread(file_name) %>%
    mutate(biosample = file_name)
  
  if(nrow(temp) > 0) {
    return(temp)
  } else {
    return(NULL)
  }
}

merged <- bind_rows(morsels) %>%
  mutate(biosample = gsub(file_dir, "", biosample)) %>%
  mutate(biosample = gsub(".csv", "", biosample)) %>%
  filter(biosample %in% names(fna))

n_distinct(merged$biosample) == length(fna) 

tibble(biosample = names(fna)) %>%
  filter(!(biosample %in% merged$biosample)) %>%
  left_join(meta %>% select(biosample, dataset)) %>%
  View()
