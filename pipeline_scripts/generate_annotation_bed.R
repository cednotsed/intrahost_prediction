rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)

df <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.csv")

parsed <- df %>%
  filter(!is.na(protein_length)) %>%
  group_by(protein_name) %>%
  summarise(start = min(nucleotide_pos), end = max(nucleotide_pos)) %>%
  mutate(ref = "MN908947.3") %>%
  mutate(protein_name = ifelse(protein_name == "NSP12a (part 1)", "NSP12_1", protein_name)) %>%
  mutate(protein_name = ifelse(protein_name == "NSP12a (part 2)", "NSP12_2", protein_name)) %>%
  select(ref, start, end, protein_name) %>%
  mutate(start = start - 1,
         end = end - 1)

fwrite(parsed, "data/metadata/wuhan-hu-1_genome_annotations_V2.bed", 
       sep = "\t",
       col.names = F,
       eol = "\n")
