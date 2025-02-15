rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

sasa_df <- fread("data/external_datasets/6vxx.filt.SASA.csv")
bio_df <- fread("results/mutation_out/all_physiochemical_stats.csv")

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv") %>%
  filter(mutation_type == "NS") %>%
  filter(protein_name == "Spike") %>%
  distinct(ref_AA, codon_number)

parsed <- sasa_df %>%
  group_by(ref_AA, codon_number, maxASA) %>%
  summarise(SASA = mean(SASA)) %>%
  rowwise() %>% 
  mutate(RSA = SASA/maxASA)

plot_df <- parsed %>%
  left_join(bio_df) %>%
  filter(protein_name == "Spike")

count_df <- plot_df %>%
  group_by(is_upper_charge) %>%
  summarise(n = n_distinct(mutation_name))

plt1 <- plot_df %>%
  ggplot(aes(x = is_upper_charge, y = delta_charge, fill = is_upper_charge)) +
  geom_boxplot() +
  geom_pwc() +
  geom_text(aes(x = is_upper_charge, y = -2.2, label = str_glue("n={n}")),
            data = count_df) +
  scale_fill_manual(values = c("grey", "indianred")) +
  labs(x = ">75 percentile of abs. charge change", 
       y = "delta charge") +
  theme_bw() +
  theme(legend.position = "none")

plt2 <- plot_df %>%
  ggplot(aes(x = is_upper_charge, y = RSA, fill = is_upper_charge)) +
  geom_boxplot() +
  geom_pwc() +
  geom_text(aes(x = is_upper_charge, y = -0.5, label = str_glue("n={n}")),
            data = count_df) +
  scale_fill_manual(values = c("grey", "indianred")) +
  labs(x = ">75 percentile of abs. charge change", 
       y = "Relative solvent accessibility") +
  theme_bw() +
  theme(legend.position = "none")

ggarrange(plt1, plt2, nrow = 1)
ggsave("results/mutation_out/charge_RSA_distribution_spike.pdf", width = 8, height = 5)
