rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggrepel)
require(ggpubr)

possible_variants <- c("G", "A", "V", "L", "I",
                       "T", "S", "M", "C", "P",
                       "F", "Y", "W", "H", "K",
                       "R", "D", "E", "N", "Q")

monthly_df <- fread("results/allele_frequency_out/observed_out/all_monthly_frequencies.csv")

mut_meta <- monthly_df %>% 
  distinct(mutation_name) %>%
  separate(mutation_name, c("protein_name", "mut"), "_", remove = F) %>%
  mutate(codon_number = parse_number(mut)) %>%
  mutate(mut = gsub("[0-9]", "", mut)) %>%
  mutate(ref_AA = substr(mut, 1, 1),
         var_AA = substr(mut, 2, 2)) %>%
  dplyr::select(mutation_name, protein_name, ref_AA, codon_number, var_AA) %>%
  mutate(region = ifelse(grepl("NSP", protein_name), "ORF1ab", protein_name))

protein_meta <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv") %>%
  filter(!is.na(protein_length)) %>%
  distinct(region, protein_length) %>%
  mutate(region = case_when(region == "ORF1ab" ~ "ORF1ab",
                            region == "S" ~ "Spike",
                            grepl("ORF", region) & region != "ORF1ab" ~ gsub("ORF", "NS", region),
                            TRUE ~ region)) %>%
  mutate(region = ifelse(region == "NS3a", "NS3", region))

orf1ab <- protein_meta %>%
  filter(region == "ORF1ab") %>%
  group_by(region) %>%
  summarise(protein_length = sum(protein_length)) %>%
  mutate(protein_length = 7097)

protein_meta <- protein_meta %>%
  filter(region != "ORF1ab") %>%
  bind_rows(orf1ab) %>%
  # Remove start and stop codons
  mutate(protein_length = protein_length - 2) %>%
  mutate(possible_mutations = 20 * protein_length)

monthly_filt <- monthly_df %>% 
  left_join(mut_meta) %>%
  filter(var_AA %in% possible_variants,
         ref_AA %in% possible_variants) %>%
  mutate(collection_month = as.Date(str_glue("{collection_month}-01")))
  
months <- deframe(monthly_filt %>% 
  distinct(collection_month))

morsels <- foreach(month_name = months) %do% {
  # month_name = months[5]
  monthly_filt %>%
    filter(collection_month <= month_name) %>%
    # group_by(region) %>%
    summarise(n = n_distinct(mutation_name)) %>%
    mutate(collection_month = month_name)
}

plot_df <- bind_rows(morsels) %>%
  left_join(protein_meta) %>%
  mutate(prop = n / possible_mutations) %>%
  mutate(region = factor(region, c("ORF1ab", "Spike", "NS3",
                                             "E", "M", "NS6", 
                                             "NS7a", "NS7b", "NS8",
                                             "N", "NS10")))

pal <- c("#DCAC51", "#CB9E88", "#B44CE7", 
         "#DEE561", "#DDAAD2", "#D5DCD9", 
         "#D1E0A2", "#897FD4", "#84B2D8", 
         "#7ADD91", "#DC6AC2", "#90EA4F", 
         "#E45B63", "#70D7D0")

plot_df %>%
  ggplot(aes(x = collection_month, y = prop, color = region)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = pal) +
  theme_bw() +
  labs(x = "Collection month", y = "Prop. mutation space")

plot_df %>%
  filter(collection_month == "2021-12-01")
  filter(prop >= 0.4) %>%
  group_by(region) %>%
  summarise(min_month = min(collection_month)) %>%
  arrange(desc(min_month))
   
  sum(protein_meta$possible_mutations)
