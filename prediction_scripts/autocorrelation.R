rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(paletteer)


morsels <- foreach(d = c("alpha", "delta", "ba1", "ba5", "xbb", "pirola")) %do% {
  stat_df <- fread(str_glue("results/prediction_out/mutation_stats/{d}.dprime_stats.csv"))
  
  rho <- cor(stat_df$prior_n, stat_df$global_n, method = "spearman")
  tibble(alias = d, rho = rho)
}

bind_rows(morsels)
