rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)

matrix_to_dnaset <- function(mat) {
  temp_mat <- apply(mat, 1, paste0, collapse = "")
  dnaset <- DNAStringSet(temp_mat, use.names = T)
  return(dnaset)
}

meta <- fread("data/external_datasets/tonkin/tonkin.metadata.csv")

file_dir <- "results/pipeline_out.tonkin/consensus_out/"
file_list <- list.files(file_dir, full.names = T)

fna <- foreach(file_name = file_list, .combine = "c") %do% {
  readDNAStringSet(file_name)
}

mat <- as.matrix(fna)

# Remove gappy sequences
prop_gaps <- apply(mat, 1,
                   function(x) {sum(x %in% c("-", "N")) / ncol(mat)})

to_remove <- prop_gaps > 0.10
n_removed <- length(prop_gaps[to_remove])

print(str_glue("Removed {n_removed} sequences"))

mat <- mat[!to_remove, ]

masked_aln <- matrix_to_dnaset(mat)

writeXStringSet(masked_aln, "data/alignments/tonkin.gap_filtered.aln")

meta %>%
  filter(sample_name %in% names(masked_aln)) %>%
  summarise(n = n())
