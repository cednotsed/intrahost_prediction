rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)

fna <- readDNAStringSet("data/alignments/reassembled.gap_filtered.pango_filtered.aln")
meta <- fread("data/metadata/all_sra_metadata.csv")
sra_meta <- fread("data/metadata/sra_metadata/all_sra_metadata.submission_dates.tsv")

parsed <- meta %>%
  filter(biosample %in% names(fna)) %>% 
  left_join(sra_meta %>% select(biosample = BioSample, submission_date = SubmissionDate)) %>% 
  mutate(submission_date = as.Date(submission_date),
         collection_date = as.Date(collection_date, "%Y-%m-%d")) %>%
  mutate(time_diff = difftime(submission_date, collection_date, units = "weeks")) %>%
  separate(geo_loc_name, c("country"), sep = "\\:\\ ", remove = F)

codon_df <- fread("results/allele_frequency_out/codon_out/missense_freq.filt.csv.gz")

codon_df <- bind_rows(morsels) %>%
  filter(id %in% names(fna)) %>%
  filter(freq < 0.5)

codon_agg <- tibble(id = names(fna)) %>%
  left_join(codon_df %>%
            group_by(id) %>%
            summarise(n = n_distinct(mutation_name))) %>%
  mutate(n = replace_na(n, 0)) 

codon_agg %>%
  summarise(median(n),
            mean(n),
            min(n), 
            max(n))

codon_agg %>%
  left_join(meta %>% select(id = biosample, dataset)) %>%
  group_by(dataset) %>%
  summarise(n = mean(n))
 
# Tonkin-Hill
tonkin_dedup <- fread("data/external_datasets/tonkin/tonkin.dedup.txt")
tonkin_fna <- readDNAStringSet("data/alignments/tonkin.gap_filtered.aln")
tonkin_fna <- tonkin_fna[tonkin_dedup$sample_name]
tonkin <- fread("results/allele_frequency_out/codon_out/tonkin/tonkin.missense_freq.filt.dedup.csv.gz")  %>%
  filter(freq < 0.5)

tibble(id = names(tonkin_fna)) %>%
  left_join(tonkin %>%
              group_by(id) %>%
              summarise(n = n_distinct(mutation_name))) %>%
  mutate(n = replace_na(n, 0)) %>%
  summarise(median(n),
            mean(n),
            min(n), 
            max(n))
  

