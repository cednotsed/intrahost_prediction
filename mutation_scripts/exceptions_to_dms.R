rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(ggpubr)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

parsed_agg <- fread("results/mutation_out/monthly_freq_aggregate.csv")

possible_variants <- c("G", "A", "V", "L", "I",
                       "T", "S", "M", "C", "P",
                       "F", "Y", "W", "H", "K",
                       "R", "D", "E", "N", "Q")

# Add scores
aa_df <- fread("data/metadata/dms/final_variant_scores.csv") %>%
  filter(target == "Wuhan-Hu-1") %>%
  mutate(mutation_name = str_glue("Spike_{mutation}"))

abx_df <- fread("data/metadata/dms/MAP_paper_antibodies_raw_data.csv") %>%
  mutate(mutation_name = str_glue("Spike_{wildtype}{site}{mutation}")) %>%
  select(mutation_name, condition, mut_escape) %>%
  group_by(mutation_name) %>%
  summarise(mean_escape = sum(mut_escape) / n())
# pivot_wider(id_cols = mutation_name, names_from = condition, values_from = mut_escape) %>%
# column_to_rownames("mutation_name")

plot_df <- parsed_agg %>%
  mutate(is_fixed = max_prop > 0.9) %>%
  inner_join(aa_df %>% select(mutation_name, delta_bind, delta_expr)) %>%
  inner_join(abx_df)

of_interest <- plot_df %>%
  filter(grepl("357|420|440|456|473", mutation_name)) %>%
  arrange(desc(mean_escape))

of_interest <- plot_df %>%
  filter(mutation_name %in% c("Spike_L455S", "Spike_F456L"))

plot_df %>%
  ggplot(aes(x = mean_escape)) +
  geom_density() +
  geom_vline(aes(xintercept = mean_escape),
             lty = "dashed",
             data = of_interest) 
  
plot_df %>% 
  arrange(desc(mean_escape)) %>%
  filter(codon_number %in% c(200, 420, 455, 456)) %>%
  filter(grepl("Spike", mutation_name)) %>% View()

plot_df %>% 
  arrange(desc(n))
  ggplot(aes(x))
