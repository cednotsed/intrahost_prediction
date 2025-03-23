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

# dataset <- "delta"

morsels <- foreach(dataset = unique(meta$alias)) %do% {
  intra_agg <- fread(str_glue("results/allele_frequency_out/codon_out/{dataset}.csv")) %>%
    filter(ref_AA %in% possible_variants) %>%
    filter(var_AA %in% possible_variants)
  
  # Prior frequency
  before_agg <- fread(str_glue("results/allele_frequency_out/observed_out/before_{dataset}.aggregate.csv")) %>%
    select(mutation_name, prior_n = global_n)
  
  # High freq mutations during dataset timeframe
  current_agg <- fread(str_glue("results/allele_frequency_out/observed_out/during_{dataset}.aggregate.csv"))
  
  # Future frequency
  future_agg <- fread(str_glue("results/allele_frequency_out/observed_out/after_{dataset}.aggregate.csv"))
  
  merged <- intra_agg %>%
    left_join(future_agg %>% dplyr::select(mutation_name, max_prop, global_n, global_prop)) %>%
    left_join(before_agg) %>%
    mutate(prior_n = replace_na(prior_n, 0)) %>%
    mutate(global_n = replace_na(global_n, 0)) %>%
    mutate(global_prop = replace_na(global_prop, 0)) %>%
    mutate(max_prop = replace_na(max_prop, 0)) %>%
    mutate(future_high = max_prop > 0.1) %>%
    mutate(future_fixed = max_prop > 0.9)
  
  # Count number of current mutations
  n_current <- merged %>%
    filter(mutation_name %in% current_agg$mutation_name) %>%
    distinct(mutation_name) %>%
    nrow()
  
  perc_current <- n_current / n_distinct(merged$mutation_name) * 100
  
  print(str_glue("{dataset}: current={n_current} ({perc_current}%)"))
  
  merged_filt <- merged %>%
    filter(!(mutation_name %in% current_agg$mutation_name)) # Remove mutations currently at high frequency
  
  merged_filt %>% 
    fwrite(str_glue("results/mutation_stats/{dataset}.stats.csv"),
           eol = "\n")
  
  table(merged_filt$future_fixed)
  table(merged_filt$future_high)
  
  # Count n_high
  n_high <- merged_filt %>%
    filter(future_high) %>%
    nrow()
  perc_high <- n_high / nrow(merged_filt) * 100
  
  tibble(dataset = dataset,
         n_high = n_high,
         perc_high = perc_high,
         n_current = n_current,
         perc_current = perc_current,
         n_total = nrow(merged),
         n_filt = nrow(merged_filt))
}

bind_rows(morsels) %>% View()
