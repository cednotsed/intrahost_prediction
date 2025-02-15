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

plot_df %>%
  ggplot(aes(x = delta_bind, y = mean_escape)) +
  geom_point() +
  geom_smooth()

cor.test(plot_df$delta_bind, plot_df$mean_escape, method = "spearman")
