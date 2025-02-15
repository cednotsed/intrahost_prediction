rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

fna <- readDNAStringSet("data/alignments/reassembled.masked.filt.aln")
meta <- fread("data/metadata/all_sra_metadata.csv")
ivar_df <- fread("results/allele_frequency_out/ivar_out/ivar_freq.all.raw.csv")
codon_df <- fread("results/allele_frequency_out/codon_out/codon_freq.all.raw.csv")

codon_filt <- codon_df %>%
  filter(id %in% names(fna)) %>%
  filter(coverage >= 100) %>%
  filter(freq > 0.03) %>%
  filter(ref_AA != "*") %>%
  filter(var_AA != "*")

ivar_filt <- ivar_df %>%
  filter(id %in% names(fna)) %>%
  filter(total_dp >= 100) %>%
  filter(alt_freq > 0.03) %>%
  filter(ref_aa != "*") %>%
  filter(alt_aa != "*") 

codon_count <- codon_filt %>%
  group_by(id) %>%
  summarise(n_codon = n_distinct(mutation_name))

ivar_count <- ivar_filt %>%
  group_by(id) %>%
  summarise(n_ivar = n_distinct(mutation_name))

plot_df <- codon_count %>%
  left_join(ivar_count)

  filter(n_ivar > 1000)
  ggplot(aes(x = n_ivar, y = n_codon)) +
  geom_point()

plot_df %>%
  summarise(median(n_codon, na.rm = T),
            median(n_ivar, na.rm = T))


