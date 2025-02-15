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

total_mutations <- 9733 * 20

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
    summarise(n = n_distinct(mutation_name)) %>%
    mutate(collection_month = month_name)
}

plot_df <- bind_rows(morsels) %>%
  mutate(prop = n / total_mutations)

plot_df %>%
  ggplot(aes(x = collection_month, y = prop)) +
  geom_line(color = "steelblue", size = 3) +
  geom_point() +
  theme_bw() +
  labs(x = "Collection month", y = "Prop. mutation space")

ggsave("results/mutation_out/prop_mutations_over_time.pdf", dpi = 600,
       height = 4, width = 4)

plot_df %>%
  filter(prop >= 0.5)
  filter(collection_month == "2022-12-01")
  
filter(prop >= 0.5) %>%
  group_by(region) %>%
  summarise(min_month = min(collection_month)) %>%
  arrange(desc(min_month))

sum(protein_meta$possible_mutations)
