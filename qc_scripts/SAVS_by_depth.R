rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)

fna <- readDNAStringSet("data/alignments/reassembled.masked.filt.aln")
meta <- fread("data/metadata/all_sra_metadata.csv")

file_dir <- "results/pipeline_out.020225/stats_out/"
file_list <- list.files(file_dir, "coverage", full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  fread(file_name) %>%
    mutate(id = file_name) %>%
    mutate(across(everything(), as.character))
}

merged <- bind_rows(morsels) %>%
  mutate(id = gsub(file_dir, "", id)) %>%
  mutate(id = gsub(".coverage.txt", "", id))

codon_df <- fread("results/allele_frequency_out/codon_out/codon_freq.all.raw.csv")

codon_filt <- codon_df %>%
  filter(id %in% names(fna)) %>%
  filter(coverage >= 100) %>%
  filter(freq > 0.03) %>%
  filter(freq < 0.5) %>%
  filter(ref_AA != "*") %>%
  filter(var_AA != "*")

codon_filt %>% 
  group_by(id) %>% 
  summarise(n = n_distinct(mutation_name)) %>%
  ungroup() %>%
  summarise(mean(n))
# ivar_df <- fread("results/allele_frequency_out/ivar_out/ivar_freq.all.raw.csv")
# ivar_filt <- ivar_df %>%
#   filter(id %in% names(fna)) %>%
#   filter(total_dp >= 100) %>%
#   filter(alt_freq > 0.03) %>%
#   filter(alt_freq < 0.5) %>%
#   filter(ref_aa != "*") %>%
#   filter(alt_aa != "*")
# ivar_filt %>%
#   group_by(id) %>%
#   summarise(n = n_distinct(mutation_name)) %>%
#   ungroup() %>%
#   summarise(mean(n))

codon_counts <- codon_filt %>%
  group_by(id) %>%
  summarise(n_codons = n_distinct(mutation_name))

plot_df <- merged %>%
  filter(id %in% names(fna)) %>% 
  mutate(numreads = as.numeric(numreads)) %>%
  left_join(codon_counts)
  
plot_df %>%
  ggplot(aes(x = log10(numreads), y = n_codons)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Log10(uniquely mapped reads)", y = "No. subconsensus SAVs")

tibble(id = names(fna)) %>%
  left_join(merged) %>% 
  filter(is.na(numreads))
