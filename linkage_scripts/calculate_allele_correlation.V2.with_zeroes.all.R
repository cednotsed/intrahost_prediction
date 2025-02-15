rm(list = ls())
setwd("/mnt/c/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)

# dataset <- "early"
meta <- fread("data/metadata/all_sra_metadata.csv")

possible_variants <- c("G", "A", "V", "L", "I",
                       "T", "S", "M", "C", "P",
                       "F", "Y", "W", "H", "K",
                       "R", "D", "E", "N", "Q")

args <- commandArgs(trailingOnly=TRUE)
dataset <- args[1]
# foreach(dataset = unique(meta$alias)) %do% {
  print(dataset)
  freq_agg <- fread(str_glue("results/allele_frequency_out/codon_out/{dataset}.csv"))
  
  freq_df <- fread(str_glue("results/allele_frequency_out/codon_out/{dataset}.raw.filt.csv")) %>%
    filter(mutation_name %in% freq_agg$mutation_name)
  
  mat <- freq_df %>%
    select(id, mutation_name, freq) %>%
    pivot_wider(id_cols = id, names_from = mutation_name, values_from = freq) %>%
    column_to_rownames("id")
  
  mat[is.na(mat)] <- 0
  
  # Get mutation counts
  col_sums <- colSums(mat > 0)
  count_df <- tibble(mutation = names(col_sums), 
                     n_biosamples = col_sums)
  to_keep <- names(col_sums)[col_sums > 5]
  
  mat_filt <- mat[, to_keep]
  
  pairs <- combn(to_keep, 2)
  
  print(ncol(pairs))
  
  morsels <- foreach(i = seq(ncol(pairs))) %do% {
    if(i %% 1000 == 0) {
      print(i)
    }
    # print(i)
    # pair <- c("N_R203K", "N_G204R")
    pair <- pairs[, i]
    mat_temp <- mat[, pair]
    mat_bool <- mat_temp > 0
    samples_to_keep <- rowSums(mat_bool) > 0
    samples_to_keep2 <- rowSums(mat_bool) == 2
    
    if(sum(samples_to_keep) >= 5) {
      corr <- cor(mat_temp[, 1], mat_temp[, 2])
      
      # Remove double zeroes
      mat_temp1 <- mat_temp[samples_to_keep, ]
      corr_no_double_zeroes <- cor(mat_temp1[, 1], mat_temp1[, 2])
      
      # Remove single zeroes
      mat_temp2 <- mat_temp[samples_to_keep2, ]
      corr_no_zeroes <- cor(mat_temp2[, 1], mat_temp2[, 2])
      
      res <- tibble(mutation1 = pair[1], mutation2 = pair[2], 
                    corr = corr, 
                    corr_no_zeroes = corr_no_zeroes,
                    corr_no_double_zeroes = corr_no_double_zeroes,
                    n_no_double_zeroes = sum(samples_to_keep),
                    n_no_zeroes = sum(samples_to_keep2),
                    mut1_min = min(mat_temp[, 1]),
                    mut1_max = max(mat_temp[, 1]),
                    mut2_min = min(mat_temp[, 2]),
                    mut2_max = max(mat_temp[, 2]),
                    mut1_cov = sd(mat_temp[, 1]) / mean(mat_temp[, 1]),
                    mut2_cov = sd(mat_temp[, 1]) / mean(mat_temp[, 2]),
                    mut1_cov_no_double = sd(mat_temp1[, 1]) / mean(mat_temp1[, 1]),
                    mut2_cov_no_double = sd(mat_temp1[, 2]) / mean(mat_temp1[, 2]),
                    mut1_cov_no_zeroes = sd(mat_temp2[, 1]) / mean(mat_temp2[, 1]),
                    mut2_cov_no_zeroes = sd(mat_temp2[, 2]) / mean(mat_temp2[, 2]))
      
      return(res)
    } else {
      return(NULL)
    }
  }
  
  merged <- bind_rows(morsels) %>%
    arrange(desc(corr)) %>%
    mutate(rsquared = corr * corr) %>%
    left_join(count_df %>% dplyr::rename(mutation1 = mutation, n_biosamples1 = n_biosamples)) %>%
    left_join(count_df %>% dplyr::rename(mutation2 = mutation, n_biosamples2 = n_biosamples)) %>%
    filter(!is.na(corr))
  
  # merged_filt <- merged %>%
  #   filter(rsquared > 0.25)
  
  fwrite(merged, str_glue("results/linkage_out/intrahost_linkage.all/{dataset}.n5.with_zeroes.csv"))
  # fwrite(merged_filt, "results/linkage_out/intrahost_linkage/linkage.early_omicron.n5.with_zeroes.rsquared25.csv")
# }
