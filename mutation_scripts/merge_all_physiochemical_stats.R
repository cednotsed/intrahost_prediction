rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

# Master df
parsed_agg <- fread("results/allele_frequency_out/observed_out/all_time.aggregate.csv")

possible_variants <- c("G", "A", "V", "L", "I",
                       "T", "S", "M", "C", "P",
                       "F", "Y", "W", "H", "K",
                       "R", "D", "E", "N", "Q")

# Add scores
aa_df <- fread("data/metadata/aa_properties.csv")

changes <- parsed_agg %>%
  distinct(ref_AA, var_AA)

# changes %>%
# filter(ref_AA == "i")
change_morsels <- foreach(i = seq(nrow(changes))) %do% {
  # i = 218
  row <- changes[i, ] %>%
    left_join(aa_df %>% select(ref_AA = amino_acid, ref_mw = mw, ref_charge = charge, ref_hydropathy = hydropathy)) %>%
    left_join(aa_df %>% select(var_AA = amino_acid, var_mw = mw, var_charge = charge, var_hydropathy = hydropathy)) %>%
    mutate(delta_mw = var_mw - ref_mw,
           delta_charge = var_charge - ref_charge,
           delta_hydropathy = var_hydropathy - ref_hydropathy)
}

change_df <- bind_rows(change_morsels) %>%
  select(ref_AA, var_AA, any_of(contains("delta"))) %>%
  mutate(abs_mw = abs(delta_mw),
         abs_charge = abs(delta_charge),
         abs_hydropathy = abs(delta_hydropathy))

# Add phenotypes
hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv") %>%
  filter(mutation_type == "NS") %>%
  distinct(protein_name)

dms_df <- fread("data/external_datasets/dms/final_variant_scores.csv") %>%
  filter(target == "Wuhan-Hu-1") %>%
  mutate(mutation_name = str_glue("Spike_{mutation}"))

abx_df <- fread("data/external_datasets/dms/MAP_paper_antibodies_raw_data.csv") %>%
  mutate(mutation_name = str_glue("Spike_{wildtype}{site}{mutation}")) %>%
  select(mutation_name, condition, mut_escape) %>%
  group_by(mutation_name) %>%
  summarise(mean_escape = sum(mut_escape) / n())

dnds_df <- fread("data/external_datasets/pond_2023.parsed.csv") %>%
  mutate(codon_name = str_glue("{protein_name}_{codon_number}"))

merged <- parsed_agg %>%
  mutate(is_fixed = max_prop > 0.9) %>%
  mutate(is_high = max_prop > 0.1) %>%
  mutate(codon_name = str_glue("{protein_name}_{codon_number}")) %>%
  left_join(change_df) %>%
  left_join(dms_df %>% select(mutation_name, delta_bind, delta_expr)) %>%
  left_join(abx_df) %>%
  left_join(dnds_df)

quantiles <- merged %>%
  summarise(upper_mw = quantile(abs_mw, 0.75),
            upper_charge = quantile(abs_charge, 0.75),
            upper_hydro = quantile(abs_hydropathy, 0.75))

plot_df <- merged %>%
  mutate(is_upper_hydro = abs_hydropathy > quantiles$upper_hydro,
         is_upper_mw = abs_mw > quantiles$upper_mw,
         is_upper_charge = abs_charge > quantiles$upper_charge) %>%
  mutate(protein_name = factor(protein_name, hookup$protein_name)) 
  # filter(!(protein_name %in% c("NS10", "NSP11")))

plot_df %>%
  fwrite("results/mutation_out/all_physiochemical_stats.csv")
