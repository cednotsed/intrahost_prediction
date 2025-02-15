rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)

# Add scores
aa_df <- fread("data/external_datasets/dms/final_variant_scores.csv") %>%
  filter(target == "Wuhan-Hu-1") %>%
  mutate(mutation_name = str_glue("Spike_{mutation}"))

abx_df <- fread("data/external_datasets/dms/MAP_paper_antibodies_raw_data.csv") %>%
  mutate(mutation_name = str_glue("Spike_{wildtype}{site}{mutation}")) %>%
  select(mutation_name, condition, mut_escape) %>%
  group_by(mutation_name) %>%
  summarise(mean_escape = sum(mut_escape, na.rm = T) / n())

aa_df %>%
  full_join(abx_df) %>%
  select(mutation_name, delta_bind, delta_expr, mean_escape) %>%
  fwrite("data/external_datasets/dms/parsed_dms_phenotypes.csv")
