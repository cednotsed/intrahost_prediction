rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

fna <- readDNAStringSet("data/alignments/tonkin.gap_filtered.aln")
df <- fread("results/allele_frequency_out/codon_out/tonkin/tonkin.missense_freq.filt.csv.gz")
meta <- fread("data/external_datasets/tonkin/tonkin.metadata.csv")
# Get intersections
duplicated <- tibble(id = names(fna)) %>%
  separate(id, c("biosample", "run"), remove = F) %>%
  group_by(biosample) %>%
  summarise(n = n_distinct(run)) %>%
  filter(n == 2)

parsed <- df %>%
  separate(id, c("biosample", "run"), remove = F) %>%
  filter(biosample %in% duplicated$biosample) %>%
  select(biosample, run, mutation_name, freq)

morsels <- foreach(bs = unique(parsed$biosample)) %do% {
  temp <- parsed %>%
    filter(biosample %in% bs)
  
  runs <- unique(temp$run)
  
  run1 <- temp %>%
    filter(run == runs[1]) %>%
    select(biosample, mutation_name, freq1 = freq)
  
  run2 <- temp %>%
    filter(run == runs[2]) %>%
    select(biosample, mutation_name, freq2 = freq)
  
  final <- run1 %>%
    full_join(run2) %>%
    mutate(freq1 = replace_na(freq1, 0),
           freq2 = replace_na(freq2, 0))
  return(final)
}

merged <- bind_rows(morsels)

# Check if intersections correlate with difference in depth
depth_meta <- meta %>%
  filter(biosample %in% unique(merged$biosample)) %>%
  group_by(biosample) %>%
  summarise(depth_ratio = max(bases) / min(bases),
            depth_diff = max(bases) - min(bases))

plot_df <- merged %>%
  mutate(both = freq1 > 0 & freq2 > 0) %>%
  group_by(biosample) %>%
  summarise(n_both = sum(both), 
            n_total = n_distinct(mutation_name)) %>%
  ungroup() %>%
  mutate(prop = n_both / n_total) %>%
  left_join(depth_meta)

cor.test(plot_df$prop, plot_df$depth_diff, method = "spearman")

plot_df %>%  
  ggplot(aes(x = depth_diff, y = prop)) +
  geom_bin2d() +
  geom_smooth()

cor(merged$freq1, merged$freq2)

# Get largest library
meta_filt <- meta %>%
  filter(sample_name %in% names(fna)) %>%
  group_by(biosample) %>%
  arrange(desc(spots)) %>%
  filter(row_number() == 1) 

meta_filt %>%
  select(sample_name) %>%
  fwrite("data/external_datasets/tonkin/tonkin.dedup.txt")

df_filt <- df %>%
  filter(id %in% meta_filt$sample_name)

df_filt %>%
  separate(id, c("biosample", "run"), "\\.") %>%
  group_by(biosample) %>%
  summarise(n = n_distinct(run)) %>%
  arrange(desc(n))

fwrite(df_filt, "results/allele_frequency_out/codon_out/tonkin/tonkin.missense_freq.filt.dedup.csv.gz")
