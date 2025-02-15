rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(paletteer)
hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

df <- fread("data/external_datasets/cov-lineages.raw.csv") %>%
  filter(variant != "")


morsels <- foreach(v = unique(df$variant)) %do% {
  temp <- df %>%
    filter(variant == v)
  
  splits <- str_split(temp$mutation_string, "\\,")[[1]]
  splits <- trimws(splits)
  
  tibble(variant = temp$variant, raw_mutations = splits)
}

merged <- bind_rows(morsels) %>% 
  separate(raw_mutations, c("protein_name", "mutation"), "\\:", remove = F) %>%
  mutate(codon_pos = extract_numeric(mutation)) %>% 
  rowwise() %>%
  mutate(mutation = gsub(unique(codon_pos), "", mutation)) %>% 
  separate(mutation, c(NA, "ref", "var"), "", remove = F) %>% 
  filter(!(protein_name %in% c("nuc", "del"))) %>%
  filter(!grepl("-", raw_mutations)) %>%
  mutate(protein_name = toupper(protein_name)) 
  
orf1ab <- merged %>%
  filter(grepl("orf1|1ab", protein_name, ignore.case = T))

non_orf1ab <- merged %>%
  filter(!(protein_name %in% orf1ab$protein_name))

n_orf1 <- (13468 - 266 + 1) / 3

orf1ab_hookup <- hookup %>% 
  filter(region == "ORF1ab") %>%
  distinct(region, codon_pos = codon_from_gene_start, new_protein_name = protein_name, codon_number)

parsed_orf1ab <- orf1ab %>%
  mutate(codon_pos = ifelse(protein_name == "ORF1b", codon_pos + n_orf1, codon_pos)) %>%
  mutate(protein_name = "ORF1ab") %>%
  left_join(orf1ab_hookup %>% dplyr::rename(protein_name = region)) %>% 
  select(variant, new_protein_name, ref, codon_number, var) %>%
  dplyr::rename(protein_name = new_protein_name, codon_pos = codon_number)

parsed_df <- bind_rows(parsed_orf1ab, non_orf1ab) %>% 
  filter(var != "*") %>%
  mutate(protein_name = ifelse(protein_name %in% c("SPIKE", "S"), "Spike", protein_name)) %>%
  mutate(protein_name = ifelse(protein_name %in% c("8", "ORF8"), "NS8", protein_name)) %>%
  mutate(protein_name = ifelse(grepl("orf3", protein_name, ignore.case = T), "NS3", protein_name)) %>%
  mutate(protein_name = ifelse(grepl("orf7a", protein_name, ignore.case = T), "NS7a", protein_name)) %>%
  mutate(protein_name = ifelse(grepl("orf6", protein_name, ignore.case = T), "NS6", protein_name)) %>%
  mutate(mutation_name = str_glue("{protein_name}_{ref}{codon_pos}{var}"))

agg_df <- fread("results/allele_frequency_out/observed_out/all_time.aggregate.csv")

final <- parsed_df %>%
  left_join(agg_df %>% select(mutation_name, global_n, max_prop))
  
mut_filt <- final %>%
  filter(!is.na(global_n)) %>%
  filter(max_prop > 0.1) %>%
  group_by(mutation_name) %>%
  summarise(n = n_distinct(variant)) %>% 
  arrange(desc(n)) %>%
  filter(n == 1)

final %>%
  select(-raw_mutations, -mutation) %>%
  filter(mutation_name %in% mut_filt$mutation_name) %>%
  fwrite("data/external_datasets/cov-constellations.parsed.csv")

final %>%
  select(-raw_mutations, -mutation) %>%
  fwrite("data/external_datasets/cov-constellations.parsed_all.csv")
  
