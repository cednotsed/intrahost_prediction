rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

meta <- fread("data/metadata/all_sra_metadata.csv")

for(d in unique(meta$alias)) {
  temp <- fread(str_glue("results/mutation_stats/{d}.stats.csv"))
  
  temp %>%
    filter(grepl("Spike", mutation_name)) %>%
    fwrite(str_glue("results/mutation_stats/{d}.spike.stats.csv"))
}

