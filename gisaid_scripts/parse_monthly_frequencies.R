rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)

file_dir <- "results/allele_frequency_out/observed_out/monthly_frequencies.V2/"
file_list <- list.files(file_dir, full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  fread(file_name) %>%
    mutate(across(everything(), as.character))
}

bind_rows(morsels) %>% 
  filter(mutation_name != "") %>%
  mutate(collection_month = str_glue("{collection_month}-01")) %>%
  fwrite("results/allele_frequency_out/observed_out/all_monthly_frequencies.csv")
