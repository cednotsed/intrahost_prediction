rm(list = ls())
setwd("/flask/scratch/matthewsp/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(doParallel)

args <- commandArgs(trailingOnly = T)
idx <- as.numeric(args[1])

parsed_agg <- fread("results/allele_frequency_out/observed_out/all_time.aggregate.csv")
parsed_filt <- parsed_agg %>%
  filter(global_n > 1000)

savs <- parsed_filt$mutation_name

aa_meta <- fread("data/metadata/gisaid/all_aa_freq.080724.csv")
aa_strings <- aa_meta$aa_substitutions

print(length(savs))

mut <- savs[idx]
print(mut)
presence <- grepl(mut, aa_strings)
  
tibble(mut = presence) %>%
    fwrite(str_glue("data/metadata/gisaid/presence_matrix.temp/{mut}.presence.csv.gz"))

