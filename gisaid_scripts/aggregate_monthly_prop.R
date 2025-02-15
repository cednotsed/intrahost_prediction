rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)

meta <- fread("data/metadata/all_sra_metadata.csv")

timeframes <- meta %>%
  group_by(dataset, alias) %>%
  summarise(start = min(collection_month),
            end = max(collection_month)) %>%
  mutate(start = start %m-% months(1),
         end = end %m+% months(1))

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
  select(mutation_name, protein_name, ref_AA, codon_number, var_AA)

monthly_parsed <- monthly_df %>%
  left_join(mut_meta) %>%
  filter(ref_AA %in% possible_variants) %>%
  filter(var_AA %in% possible_variants)

# Get future estimates
foreach(d = unique(timeframes$alias)) %do% {
  timeframe <- timeframes %>% filter(alias == d)
  
  future <- monthly_parsed %>%
    filter(collection_month > timeframe$end)
  
  global_total <- deframe(future %>% 
    distinct(collection_month, n_total) %>%
    summarise(sum(n_total)))
  
  future_agg <- future %>%
    group_by(mutation_name) %>%
    summarise(global_n = sum(n_present),
              median_prop = median(prop),
              global_prop = sum(n_present) / global_total,
              max_prop = max(prop))
  
  current_agg <- monthly_parsed %>% 
    filter(collection_month > timeframe$start & 
             collection_month < timeframe$end) %>%
    group_by(mutation_name) %>%
    summarise(max_prop = max(prop),
              global_n = sum(n_present)) %>%
    ungroup() %>%
    filter(max_prop > 0.1)
  
  # For prior linkage calculations
  if(d != "early") {
    before <- monthly_parsed %>%
      filter(collection_month < timeframe$start) %>%
      group_by(mutation_name) %>%
      summarise(max_prop = max(prop),
                global_n = sum(n_present))
    fwrite(before, str_glue("results/allele_frequency_out/observed_out/before_{d}.aggregate.csv"))
  } else {
    before <- monthly_parsed %>%
      filter(collection_month < timeframe$end) %>%
      group_by(mutation_name) %>%
      summarise(max_prop = max(prop),
                global_n = sum(n_present))
    fwrite(before, str_glue("results/allele_frequency_out/observed_out/before_{d}.aggregate.csv"))
  }
  
  fwrite(future, str_glue("results/allele_frequency_out/observed_out/after_{d}.monthly.csv"))
  fwrite(future_agg, str_glue("results/allele_frequency_out/observed_out/after_{d}.aggregate.csv"))
  fwrite(current_agg, str_glue("results/allele_frequency_out/observed_out/during_{d}.aggregate.csv"))
}

