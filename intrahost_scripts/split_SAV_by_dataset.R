rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

codon_df <- fread("results/allele_frequency_out/codon_out/missense_freq.filt.csv.gz")

# Add physiochemistry and DMS scores
aa_df <- fread("data/metadata/aa_properties.blosum62.csv")
dms_df <- fread("data/external_datasets/dms/parsed_dms_phenotypes.csv")

# Split datasets
fna <- readDNAStringSet("data/alignments/reassembled.gap_filtered.pango_filtered.aln")
meta <- fread("data/metadata/all_sra_metadata.csv") %>%
  filter(biosample %in% names(fna)) %>%
  mutate(save_prefix = tolower(gsub("\\.", "", dataset)))

foreach(d = unique(meta$dataset)) %do% {
  temp <- meta %>%
    filter(dataset == d)
  
  save_prefix <- unique(temp$save_prefix)
  
  dataset_temp <- codon_df %>%
    filter(id %in% temp$biosample)
  
  dataset_agg <- dataset_temp %>%
    group_by(protein_name, mutation_name, ref_AA, codon_number, var_AA) %>%
    summarise(n = n_distinct(id),
              n_filt = n_distinct(id),
              median_freq = median(freq),
              max_freq = max(freq),
              max_codon_variants = max(n_codons),
              median_codon_variants = median(n_codons)) %>%
    left_join(aa_df) %>%
    left_join(dms_df) 
  
  fwrite(dataset_temp, str_glue("results/allele_frequency_out/codon_out/{save_prefix}.raw.filt.csv"))
  fwrite(dataset_agg, str_glue("results/allele_frequency_out/codon_out/{save_prefix}.csv"))
}
