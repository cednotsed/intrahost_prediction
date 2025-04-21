rm(list = ls())
setwd("/flask/scratch/matthewsp/intrahost_prediction")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(doParallel)

aa_meta <- fread("data/metadata/gisaid/all_aa_freq.080724.csv") %>%
  select(collection_month)

file_dir <- "data/metadata/gisaid/presence_matrix.temp/"
file_list <- list.files(file_dir, full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  # file_name = file_list[1]
  mut <- gsub(file_dir, "", file_name)
  mut <- gsub(".presence.csv.gz", "", mut)
  mut <- gsub("\\/", "", mut)
  
  print(mut)
  temp <- fread(file_name)
  colnames(temp) <- mut
  gc()
  return(temp)
}

mat <- bind_cols(aa_meta, morsels)

fwrite(mat, "data/metadata/gisaid/presence_absence_matrix.csv")
mat[1:5, 1:5]


