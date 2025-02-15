rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)

fna <- readDNAStringSet("data/alignments/reassembled.masked.filt.aln")

file_dir <- "results/pipeline_out.020225/stats_out/"
file_list <- list.files(file_dir, "coverage", full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  id <- gsub(file_dir, "", file_name)
  id <- gsub("\\/", "", id)
  
  temp <- fread(file_name)
  
  tibble(id = id, n_rows = nrow(temp))
}

merged <- bind_rows(morsels) %>%
  mutate(id = gsub(".coverage.txt", "", id))

merged_filt <- merged %>%
  filter(n_rows > 0)

head(merged)

print(names(fna)[!(names(fna) %in% merged_filt$id)])
# tibble(names(fna)[!(names(fna) %in% merged_filt$id)]) %>%
#   fwrite("test.txt")
# 
# length(names(fna)[!(names(fna) %in% merged_filt$id)])
