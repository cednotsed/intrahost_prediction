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

meta <- fread("data/metadata/all_sra_metadata.csv")

file_dir <- "results/pipeline_out.020225/consensus_out/"
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

# # Mask gappy sites
# prop_site_gaps <- apply(mat, 2,
#                         function(x) {sum(x %in% c("-", "N")) / nrow(mat)})
# 
# site_to_mask <- prop_site_gaps > 0.10
# 
# mat[, site_to_mask] <- "N"

masked_aln <- matrix_to_dnaset(mat)

writeXStringSet(masked_aln, "data/alignments/reassembled.gap_filtered.aln")

meta %>%
  filter(biosample %in% names(masked_aln)) %>%
  group_by(dataset) %>%
  summarise(n = n())

# Remove after tempest analysis
# masked_filt <- masked_aln[!(names(masked_aln) %in% c("SAMN29663258","SAMN40189392", "SAMN26809266",
#                                                      "SAMN26809252", "SAMN39462117", "SAMN26197536"))]
# masked_filt
# writeXStringSet(masked_aln, "data/alignments/reassembled.all.masked.tempest1.aln")

# # snp_df <- bind_rows(fread("results/pipeline_out.transition/snp_counts.transition.csv"),
# #                     fread("results/pipeline_out/snp_counts.csv")) %>%
# #   filter(genome_name != "MN908947.3") %>%
# #   filter(genome_name %in% meta$biosample) %>%
# #   left_join(meta %>% select(genome_name = biosample, collection_date, collection_month, dataset)) %>%
# #   mutate(collection_date = as.Date(collection_date, "%Y-%m-%d")) %>%
# #   # filter(!(dataset == "early" & snps > 10)) %>%
# #   filter(dataset %in% c("early", "early_delta", "early_omicron"))
# 
# snp_df %>%
#   ggplot(aes(x = collection_date, y = snps)) +
#   geom_point() +
#   geom_smooth(method = "lm")
#
# fna <- fna[names(fna) %in% snp_df$genome_name]
