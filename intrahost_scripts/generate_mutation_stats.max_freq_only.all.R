rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

meta <- fread("data/metadata/all_sra_metadata.csv")
set.seed(66)

for(d in unique(meta$alias)) {
  temp <- fread(str_glue("results/mutation_stats/{d}.stats.csv"))
  colnames(temp)   
  temp$delta_charge <- sample(temp$delta_charge)
  temp$abs_charge <- sample(temp$abs_charge)
  temp$delta_mw <- sample(temp$delta_mw)
  temp$abs_mw <- sample(temp$abs_mw)
  temp$delta_hydropathy <- sample(temp$delta_hydropathy)
  temp$abs_hydropathy <- sample(temp$abs_hydropathy)
  temp$delta_bind <- sample(temp$delta_bind)
  temp$delta_expr <- sample(temp$delta_expr)
  temp$mean_escape <- sample(temp$mean_escape)
  temp$blosum62_score <- sample(temp$blosum62_score)
  
  temp %>%
    fwrite(str_glue("results/mutation_stats/{d}.max_freq_only.stats.csv"))
}

