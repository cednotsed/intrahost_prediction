rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
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

d <- "delta"

monthly_df <- fread("results/allele_frequency_out/observed_out/all_monthly_frequencies.csv") %>%
  mutate(collection_date = as.Date(str_glue("{collection_month}-01")))

time_morsels <- foreach(d = unique(meta$alias)) %do% {
  timeframes <- meta %>%
    group_by(dataset, alias) %>%
    summarise(start = min(collection_month),
              end = max(collection_month)) %>%
    mutate(start = start %m-% months(1),
           end = end %m+% months(1)) %>%
    mutate(end_month = as.numeric(format(end, "%m")),
           end_year = as.numeric(format(end, "%Y")))
  
  timeframe <- timeframes %>%
    filter(alias == d)
  
  mutation_stats <- fread(str_glue("results/modelling_out/mutation_stats/{d}.stats.csv")) %>%
    filter(future_high)
  
  monthly_filt <- monthly_df %>%
    filter(collection_date > timeframe$end)
  
  morsels <- foreach(mut = unique(mutation_stats$mutation_name)) %do% {
    monthly_filt %>%
      filter(mutation_name == mut) %>%
      arrange(desc(prop)) %>%
      head(1)
  }
  
  bind_rows(morsels) %>%
    select(mutation_name, collection_month) %>%
    separate(collection_month, c("year", "month"), "-") %>%
    mutate(month_diff = (12 * as.numeric(year) + as.numeric(month)) - (12 * timeframe$end_year + timeframe$end_month)) %>%
    mutate(alias = d)
}

merged <- bind_rows(time_morsels)

merged %>%
  group_by(alias) %>%
  summarise(median_lag = median(month_diff)) %>%
  arrange(desc(median_lag))
