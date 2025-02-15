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

timepoints <- list(c("early", "delta"),
                   c("delta", "xbb"),
                   c("early", "ba1"),
                   c("ba5", "end"),
                   c("alpha", "ba1"),
                   c("ba1", "xbb"),
                   c("xbb", "end"),
                   c("delta", "ba5"),
                   c("ba1", "ba5"),
                   c("xbb", "pirola"),
                   c("ba5", "xbb"),
                   c("pirola", "end"),
                   c("ba1", "end"),
                   c("alpha", "delta"),
                   c("delta", "ba1"))

foreach(timepoint = timepoints) %do% {
  start_point <- timepoint[1]
  end_point <- timepoint[2]
  
  print(str_glue("{start_point}_to_{end_point}"))
  
  # Get future estimates
  start_timeframe <- timeframes %>% filter(alias == start_point)
  end_timeframe <- timeframes %>% filter(alias == end_point)
  
  if(end_point != "end") {
    future <- monthly_parsed %>%
      filter(collection_month > start_timeframe$end,
             collection_month < end_timeframe$start)
  } else {
    future <- monthly_parsed %>%
      filter(collection_month > start_timeframe$end)
  }
  
  global_total <- deframe(future %>% 
                            distinct(collection_month, n_total) %>%
                            summarise(sum(n_total)))
  
  time_df <- future %>%
    summarise(start = min(collection_month),
              end = max(collection_month))
  
  future_agg <- future %>%
    group_by(mutation_name) %>%
    summarise(global_n = sum(n_present),
              median_prop = median(prop),
              global_prop = sum(n_present) / global_total,
              max_prop = max(prop)) %>%
    bind_cols(time_df)
  
  fwrite(future_agg, str_glue("results/allele_frequency_out/observed_out/{start_point}_to_{end_point}.aggregate.csv"))
}
