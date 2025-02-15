rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)

file_dir <- "results/linkage_out/observed_Dprime.temp/"

# dataset <- "pirola"

foreach(dataset = c("alpha", "delta", "ba1", "ba5", "xbb", "pirola")) %do% {
  file_list <- list.files(file_dir, full.names = T)
  file_list <- file_list[grepl(dataset, file_list)]
  
  print(length(file_list))
  
  morsels <- foreach(file_name = file_list) %do% {
    fread(file_name)  
  }
  
  bind_rows(morsels) %>%
    fwrite(str_glue("results/linkage_out/observed_linkage/observed_Dprime.before_{dataset}.gt1000.csv.gz"))
}
