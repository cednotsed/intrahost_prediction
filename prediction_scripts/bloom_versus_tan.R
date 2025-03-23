rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

bloom_df <- fread("data/external_datasets/aamut_fitness_all.bloom2023.csv") %>%
  mutate()

gisaid_df <- fread("results/allele_frequency_out/observed_out/all_time.aggregate.csv")
gisaid_df %>% distinct(protein_name)
parsed <- bloom_df %>%
  # distinct(gene) %>%
  mutate(protein_name = gsub("nsp", "NSP", gene)) %>%
  mutate(protein_name = gsub("\\ \\(Mpro\\)|\\ \\(RdRp\\)", "", protein_name)) %>%
  mutate(protein_name = ifelse(protein_name == "S", "Spike", protein_name)) %>%
  mutate(protein_name = gsub("ORF", "NS", protein_name)) %>%
  mutate(protein_name = gsub("NS3a", "NS3", protein_name)) %>% 
  mutate(mutation_name = str_glue("{protein_name}_{aa_mutation}")) %>%
  select(protein_name, mutation_name, delta_fitness) %>% 
  left_join(gisaid_df %>% select(mutation_name, global_n)) %>% 
  filter(protein_name != "NS1ab") %>%
  filter(!is.na(global_n))
  
parsed
# filter(is.na(global_n)) %>% 
  # group_by(protein_name) %>%
  # summarise(n = n()) 

parsed %>%
  ggplot(aes(x = delta_fitness, y = log10(global_n + 1))) +
  geom_point() +
  labs(x = "Fitness effect (Bloom and Neher)", y = "Log10(Mutational fitness + 1) (Tan et al.)")

cor.test(parsed$delta_fitness, log10(parsed$global_n + 1), method = "spearman") 

bloom_df %>% 
  filter(gene == "ORF1ab")

gisaid_df %>%
  distinct(protein_name)

gisaid_df %>% View()
