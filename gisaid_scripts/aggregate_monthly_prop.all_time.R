rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)

possible_variants <- c("G", "A", "V", "L", "I",
                       "T", "S", "M", "C", "P",
                       "F", "Y", "W", "H", "K",
                       "R", "D", "E", "N", "Q")

monthly_df <- fread("results/allele_frequency_out/observed_out/all_monthly_frequencies.csv") %>%
  filter(mutation_name != "") %>%
  filter(!grepl("del|ins|stop", mutation_name))

mut_meta <- monthly_df %>% 
  distinct(mutation_name) %>%
  separate(mutation_name, c("protein_name", "mut"), "_", remove = F) %>%
  mutate(codon_number = parse_number(mut)) %>%
  mutate(mut = gsub("[0-9]", "", mut)) %>%
  mutate(ref_AA = substr(mut, 1, 1),
         var_AA = substr(mut, 2, 2)) %>%
  dplyr::select(mutation_name, protein_name, ref_AA, codon_number, var_AA) %>%
  mutate(region = ifelse(grepl("NSP", protein_name), "ORF1ab", protein_name))

monthly_parsed <- monthly_df %>%
  left_join(mut_meta) %>%
  filter(ref_AA %in% possible_variants) %>%
  filter(var_AA %in% possible_variants)

global_total <- deframe(monthly_parsed %>% 
                          distinct(collection_month, n_total) %>%
                          summarise(sum(n_total)))
  
parsed_agg <- monthly_parsed %>%
  group_by(mutation_name, protein_name, region, ref_AA, codon_number, var_AA) %>%
  summarise(global_n = sum(n_present),
            global_prop = sum(n_present) / global_total,
            max_prop = max(prop))
  
fwrite(parsed_agg, str_glue("results/allele_frequency_out/observed_out/all_time.aggregate.csv"))
