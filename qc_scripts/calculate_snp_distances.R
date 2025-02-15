rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)

mat <- as.matrix(readDNAStringSet("data/alignments/reassembled.gap_filtered.aln"))

# Compute distances
ref_vec <- as.matrix(readDNAStringSet("data/genomes/MN908947.3.fna"))[1, ]

# Number of segregating sites
seg_sites <- apply(mat, 1, function(x) {sum(x != ref_vec)})

# No. of Ns 
N_sites <- apply(mat, 1, function(x) {sum(x == "N")})

# SNP count = seg. sites - N sites
snp_df <- tibble(biosample = names(seg_sites),
       n_seg = seg_sites) %>%
  left_join(tibble(biosample = names(N_sites), 
                   n_Ns = N_sites)) %>%
  mutate(n_SNPs = n_seg - n_Ns)

snp_df %>%
  fwrite("results/pipeline_out.020225/snp_counts.csv")

plot_df <- snp_df %>%
  left_join(meta) %>%
  mutate(collection_date = as.Date(collection_date, "%Y-%m-%d"))

plot_df %>%  
  ggplot(aes(x = as.Date(collection_date, "%Y-%m-%d"), y = n_SNPs)) +
  geom_point() +
  geom_smooth(method = "lm")

ggsave("results/qc_out/SNPs_v_collection_date.pdf", dpi = 600, width = 4, height = 4)
