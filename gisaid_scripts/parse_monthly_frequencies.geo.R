rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)

for(region in c("asia", "africa", "north_america", 
                "south_america", "europe", "US_only", 
                "non_US", "oceania")) {
  file_dir <- str_glue("results/allele_frequency_out/observed_out/geographical_splits/monthly_frequencies.{region}/")
  file_list <- list.files(file_dir, full.names = T)

  morsels <- foreach(file_name = file_list) %do% {
    fread(file_name) %>%
      mutate(across(everything(), as.character))
  }
  
  bind_rows(morsels) %>% 
    filter(mutation_name != "") %>%
    mutate(collection_month = str_glue("{collection_month}-01")) %>%
    fwrite(str_glue("results/allele_frequency_out/observed_out/geographical_splits/{region}.all_monthly_frequencies.csv"))
}
