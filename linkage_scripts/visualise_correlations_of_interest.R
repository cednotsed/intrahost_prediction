rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

freq_df <- fread("results/allele_frequency_out/codon_out/early.raw.filt.csv") %>%
  select(id, mutation_name, freq)

mat <- freq_df %>%
  pivot_wider(id_cols = id, names_from = mutation_name, values_from = freq) %>%
  column_to_rownames("id") %>%
  mutate(across(everything(), ~replace_na(., 0)))

mut1 <- "Spike_D614G"
mut2 <- "NSP12_P323L"

mut1 <- "N_R203K"
mut2 <- "N_G204R"

mut1 <- "NSP4_K35R"
mut2 <- "NSP4_S34F"

corr <- cor.test(mat[[mut1]], mat[[mut2]])
r <- signif(corr$estimate, 2)
p <- signif(corr$p.value, 2)

mat %>%
  ggplot(aes(x = get(mut1), y = get(mut2))) +
  geom_point(fill = "darkviolet", color = "black", alpha = 0.5,
             size = 3,
             pch = 21) +
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(slope = 1, intercept = 0, lty = "dashed") +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        axis.title = element_text(face = "bold")) +
  labs(x = str_glue("Intrahost freq. ({mut1})"), y = str_glue("Intrahost freq. ({mut2})"), title = str_glue("r = {r}, p = {p}"))

ggsave(str_glue("results/linkage_out/{mut1}_{mut2}.pdf"), dpi = 600, width = 4, height = 4)
