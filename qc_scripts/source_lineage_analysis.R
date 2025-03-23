rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)

fna <- readDNAStringSet("data/alignments/reassembled.gap_filtered.pango_filtered.aln")
meta <- fread("data/metadata/all_sra_metadata.csv")
const_df <- fread("data/external_datasets/cov-constellations.parsed_all.csv") %>%
  dplyr::rename(source_variant = variant)
sra_meta <- fread("data/metadata/sra_metadata/all_sra_metadata.submission_dates.tsv")
 
parsed <- meta %>%
  filter(biosample %in% names(fna)) %>% 
  left_join(sra_meta %>% select(biosample = BioSample, submission_date = SubmissionDate)) %>% 
  mutate(submission_date = as.Date(submission_date),
         collection_date = as.Date(collection_date, "%Y-%m-%d")) %>%
  mutate(time_diff = difftime(submission_date, collection_date, units = "weeks")) %>%
  separate(geo_loc_name, c("country"), sep = "\\:\\ ", remove = F)

codon_df <- fread("results/allele_frequency_out/codon_out/missense_freq.filt.csv.gz") %>%
  filter(freq < 0.5) %>%
  left_join(const_df)

mut_count <- const_df %>%
  group_by(source_variant) %>%
  summarise(n_total = n_distinct(mutation_name))

const_counts <- const_df %>%
  group_by(source_variant) %>%
  summarise(total_const_n = n_distinct(mutation_name)) %>%
  arrange(desc(total_const_n))

# lineage_counts <- codon_df %>%
#   group_by(id) %>%
#   summarise(n_lineages = n_distinct(source_variant, na.rm = T))
# 
# lineage_counts <- tibble(id = names(fna)) %>%
#   left_join(lineage_counts) %>%
#   mutate(n_lineages = replace_na(n_lineages, 0))
# 
# single_lineage <- lineage_counts %>%
#   filter(n_lineages == 1)

plot_df <- codon_df %>% 
  # filter(id %in% single_lineage$id) %>%
  group_by(id) %>%
  summarise(n_muts = n_distinct(mutation_name),
            n_var = sum(!is.na(source_variant)),
            n_sources = n_distinct(source_variant, na.rm = T)) %>%
  ungroup()

plot_df %>%
  summarise(median_n = median(n_var),
            min_n = min(n_var),
            max_n = max(n_var))

plot_df %>%
  nrow()
  filter(n_var > 1) %>%
  nrow()

plt1 <- plot_df %>%
  ggplot(aes(x = n_var)) +
  geom_histogram(color = "black") +
  geom_vline(xintercept = stat_df$median_n,
             lty = "dashed",
             color = "red") +
  theme_bw() +
  labs(x = "No. lineage-defining subconsensus mutations", y = "No. samples")

plt2 <- plot_df %>%
  ggplot(aes(x = n_var, y = n_sources)) +
  geom_point() +
  theme_bw() +
  labs(x = "No. lineage-defining subconsensus SAVs", y = "No. source variants")


# 
# codon_df %>%
#   group_by(id) %>%
#   filter(!is.na(source_variant)) %>%
#   summarise(n_var = n_distinct(mutation_name),
#             n_sources = n_distinct(source_variant, na.rm = T)) %>%
#   arrange(desc(n_var))
# 
#   group_by(id, source_variant) %>%
#   summarise(n_var = n_distinct(mutation_name)) %>%
#   filter(!is.na(source_variant)) %>%
#   arrange(desc(n_var)) %>%
#   left_join(const_counts) %>%
#   mutate(prop = n_var / total_const_n) %>%
#   ggplot(aes(x = prop)) + 
#   geom_histogram()
# 
# multi_mutants <- plot_df %>%
#   filter(n_var > 1) %>%
#   arrange(desc(n_var))
# 
# plot_df %>%
#   filter(id %in% multi_mutants$id) %>%
#   arrange(desc(prop))
#   filter(n_sources == 1) %>% 
#   left_join(meta %>% select(id = biosample, dataset)) %>% 
#   ggplot(aes(x = n_sources)) +
#   geom_histogram()
# 
# plt3 <- plot_df2 %>%
#   ggplot(aes(x = prop)) +
#   geom_histogram(color = "black") +
#   geom_vline(xintercept = stat_df$median_prop,
#              lty = "dashed",
#              color = "red") +
#   theme_bw() +
#   labs(x = "Prop. lineage-defining subconsensus mutations", y = "No. samples") 
#   
# plot_df2 %>%
#   summarise(median(prop) * 100,
#             min(prop) * 100,
#             max(prop) * 100)

ggarrange(plt1, plt2)

cor.test(plot_df$n_var, plot_df$n_sources)
ggsave("results/qc_out/source_lineage_analysis.pdf", dpi = 600, width = 7, height = 4)
