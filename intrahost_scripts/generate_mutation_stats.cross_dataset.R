rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

possible_variants <- c("G", "A", "V", "L", "I",
                       "T", "S", "M", "C", "P",
                       "F", "Y", "W", "H", "K",
                       "R", "D", "E", "N", "Q")

meta <- fread("data/metadata/all_sra_metadata.csv")

prefixes <- c("early_to_delta", "delta_to_xbb", "early_to_ba1",
              "ba5_to_end", "alpha_to_ba1", "ba1_to_xbb",
              "xbb_to_end", "delta_to_ba5", "ba5_to_xbb", 
              "pirola_to_end", "xbb_to_pirola", "ba1_to_ba5", 
              "ba1_to_end", "alpha_to_delta", "delta_to_ba1")

foreach(prefix = prefixes) %do% {
  dataset <- str_split(prefix, "_")[[1]][1]
    
  intra_agg <- fread(str_glue("results/allele_frequency_out/codon_out/{dataset}.csv")) %>%
    filter(ref_AA %in% possible_variants) %>%
    filter(var_AA %in% possible_variants)
  
  # Prior frequency
  before_agg <- fread(str_glue("results/allele_frequency_out/observed_out/before_{dataset}.aggregate.csv")) %>%
    select(mutation_name, prior_n = global_n)
  
  # High freq mutations during dataset timeframe
  current_agg <- fread(str_glue("results/allele_frequency_out/observed_out/during_{dataset}.aggregate.csv"))
  
  # Future frequency
  future_agg <- fread(str_glue("results/allele_frequency_out/observed_out/{prefix}.aggregate.csv"))
  
  merged <- intra_agg %>%
    left_join(future_agg %>% dplyr::select(mutation_name, max_prop, global_n, global_prop)) %>%
    left_join(before_agg) %>%
    mutate(prior_n = replace_na(prior_n, 0)) %>%
    mutate(global_n = replace_na(global_n, 0)) %>%
    mutate(global_prop = replace_na(global_prop, 0)) %>%
    mutate(max_prop = replace_na(max_prop, 0)) %>%
    mutate(future_high = max_prop > 0.1) %>%
    mutate(future_fixed = max_prop > 0.9) %>% 
    filter(!(mutation_name %in% current_agg$mutation_name))
  
  merged %>% fwrite(str_glue("results/mutation_stats/{prefix}.stats.csv"))
  
  # table(merged$future_fixed)
  table(merged$future_high)
}
