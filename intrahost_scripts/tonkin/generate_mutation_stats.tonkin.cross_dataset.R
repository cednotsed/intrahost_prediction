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

meta <- fread("data/external_datasets/tonkin/tonkin.metadata.csv")

# Add physiochemistry and DMS scores
aa_df <- fread("data/metadata/aa_properties.blosum62.csv")
dms_df <- fread("data/external_datasets/dms/parsed_dms_phenotypes.csv")

intra_agg <- fread(str_glue("results/allele_frequency_out/codon_out/tonkin/tonkin.missense_freq.filt.dedup.csv.gz")) %>%
  filter(ref_AA %in% possible_variants) %>%
  filter(var_AA %in% possible_variants) %>%
  group_by(protein_name, mutation_name, ref_AA, codon_number, var_AA) %>%
  summarise(n = n_distinct(id),
            n_filt = n_distinct(id),
            median_freq = median(freq),
            max_freq = max(freq),
            max_codon_variants = max(n_codons),
            median_codon_variants = median(n_codons)) %>%
  left_join(aa_df) %>%
  left_join(dms_df)

# Prior frequency
before_agg <- fread(str_glue("results/allele_frequency_out/observed_out/before_tonkin.aggregate.csv")) %>%
  select(mutation_name, prior_n = global_n)

# High freq mutations during dataset timeframe
current_agg <- fread(str_glue("results/allele_frequency_out/observed_out/during_tonkin.aggregate.csv"))

# Future frequency
future_agg <- fread(str_glue("results/allele_frequency_out/observed_out/tonkin_to_delta.aggregate.csv"))

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
  fwrite(str_glue("results/mutation_stats/tonkin_to_delta.stats.csv"),
         eol = "\n")

table(merged$future_fixed)
table(merged$future_high)

