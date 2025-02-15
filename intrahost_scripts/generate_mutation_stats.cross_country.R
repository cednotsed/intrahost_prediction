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

dataset <- "early_to_delta"
region <- "non_US"

# dataset <- "ba1_to_end"
# region <- "US_only"

start_timeframe <- str_split(dataset, "\\_")[[1]][1]
end_timeframe <- str_split(dataset, "\\_")[[1]][3]

intra_agg <- fread(str_glue("results/allele_frequency_out/codon_out/{start_timeframe}.{region}.csv")) %>%
  filter(ref_AA %in% possible_variants) %>%
  filter(var_AA %in% possible_variants)

# Prior frequency
before_agg <- fread(str_glue("results/allele_frequency_out/observed_out/geographical_splits/before_{start_timeframe}.{region}.aggregate.csv")) %>%
  select(mutation_name, prior_n = global_n)

# High freq mutations during dataset timeframe
current_agg <- fread(str_glue("results/allele_frequency_out/observed_out/geographical_splits/during_{start_timeframe}.{region}.aggregate.csv"))

# Future frequency
future_agg <- fread(str_glue("results/allele_frequency_out/observed_out/geographical_splits/{start_timeframe}_to_{end_timeframe}.{region}.aggregate.csv"))

merged <- intra_agg %>%
  left_join(future_agg %>% dplyr::select(mutation_name, max_prop, global_n, global_prop)) %>%
  left_join(before_agg) %>%
  mutate(prior_n = replace_na(prior_n, 0)) %>%
  mutate(global_n = replace_na(global_n, 0)) %>%
  mutate(global_prop = replace_na(global_prop, 0)) %>%
  mutate(max_prop = replace_na(max_prop, 0)) %>%
  mutate(future_high = max_prop > 0.1) %>%
  mutate(future_fixed = max_prop > 0.9) %>% 
  filter(!(mutation_name %in% current_agg$mutation_name)) # Remove mutations currently at high frequency

merged %>% 
  fwrite(str_glue("results/mutation_stats/{dataset}.{region}.stats.csv"),
         eol = "\n")

table(merged$future_fixed)
table(merged$future_high)
