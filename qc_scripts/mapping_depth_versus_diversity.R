rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

fna <- readDNAStringSet("data/alignments/reassembled.gap_filtered.pango_filtered.aln")

file_dir <- "results/pipeline_out.020225/stats_out/"
file_list <- list.files(file_dir, full.names = T)
file_list <- file_list[grepl("coverage", file_list)]

morsels <- foreach(file_name = file_list) %do% {
  id <- gsub(".coverage.txt", "", file_name)
  id <- gsub(file_dir, "", id)
  
  fread(file_name) %>%
    mutate(id = id) %>%
    mutate(across(everything(), as.character))
}

codon_df <- fread("results/allele_frequency_out/codon_out/missense_freq.filt.csv.gz")

codon_parsed <- codon_df %>%
  filter(freq < 0.5) %>%
  group_by(id) %>%
  summarise(n = n())

merged <- bind_rows(morsels) %>%
  filter(id %in% names(fna)) %>%
  select(id, numreads, meandepth, meanbaseq) %>%
  mutate(numreads = as.numeric(numreads),
         meandepth = as.numeric(meandepth),
         meanbaseq = as.numeric(meanbaseq)) %>%
  left_join(codon_parsed) %>%
  mutate(n = replace_na(n, 0))

# merged %>% 
#   ggplot(aes(x = n_alleles)) +
#   geom_histogram()
merged %>%
  ggplot(aes(x = log10(numreads), y = n)) +
  geom_point() +
  geom_smooth()

merged %>%
  ggplot(aes(x = log10(meandepth), y = n)) +
  geom_point() +
  geom_smooth()

merged %>%
  ggplot(aes(x = meanbaseq, y = n)) +
  geom_point() +
  geom_smooth(method = "lm")

cor.test(merged$numreads, merged$n, method = "spearman")
cor.test(merged$meandepth, merged$n, method = "spearman")
cor.test(merged$meanbaseq, merged$n, method = "spearman")

merged %>%
  summarise(median_depth = median(meandepth),
            median_baseq = median(meanbaseq),
            median_reads = median(numreads))

merged %>%

merged %>%
  summarise(prop =  sum(numreads > 50000) / n())
  filter(numreads > 50000)
  