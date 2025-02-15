rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

bloom_df <- fread("data/external_datasets/aamut_fitness_all.bloom2023.csv")
gisaid_df <- fread("results/allele_frequency_out/observed_out/all_time.aggregate.csv")
hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv") %>%
  filter(mutation_type == "NS") %>%
  distinct(protein_name, ref_AA, codon_number)
phy_df <- fread("data/metadata/aa_properties.blosum62.csv")

parsed <- bloom_df %>%
  # distinct(gene) %>%
  mutate(protein_name = gsub("nsp", "NSP", gene)) %>%
  mutate(protein_name = gsub("\\ \\(Mpro\\)|\\ \\(RdRp\\)", "", protein_name)) %>%
  mutate(protein_name = ifelse(protein_name == "S", "Spike", protein_name)) %>%
  mutate(protein_name = gsub("ORF", "NS", protein_name)) %>%
  mutate(protein_name = gsub("NS3a", "NS3", protein_name)) %>% 
  mutate(mutation_name = str_glue("{protein_name}_{aa_mutation}")) %>%
  select(protein_name, aa_site, mutation_name, delta_fitness, mutant_aa) %>% 
  left_join(gisaid_df %>% select(mutation_name, global_n)) %>% 
  left_join(hookup %>% select(ref_AA, protein_name, aa_site = codon_number)) %>%
  filter(protein_name != "NS1ab") %>%
  filter(!is.na(global_n)) %>%
  left_join(phy_df %>% dplyr::rename(mutant_aa = var_AA))
  
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

plt1 <- parsed %>%
  ggplot(aes(x = delta_fitness)) +
  geom_histogram(color = "black") +
  theme_classic() +
  labs(x = "Fitness effects (Bloom and Neher)", y= "No. SAVs")

plt2 <- parsed %>%
  ggplot(aes(x = factor(blosum62_score), y = delta_fitness)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "BLOSUM62 score", y = "Fitness effects (Bloom and Neher)")

ggarrange(plt1, plt2, nrow = 1)
ggsave("results/qc_out/bloom_fitness_histogram.pdf", width = 6, height = 4)
