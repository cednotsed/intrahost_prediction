rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

fna <- readDNAStringSet("data/alignments/reassembled.masked.filt.aln")
sub_df <- fread("data/metadata/all_sra_metadata.submission_dates.tsv")
const_df <- fread("data/external_datasets/cov-constellations.parsed.csv") %>%
  dplyr::rename(source_variant = variant)

meta <- fread("data/metadata/all_sra_metadata.csv") %>%
  filter(biosample %in% names(fna)) %>%
  left_join(sub_df %>% select(biosample = BioSample, submission_date = SubmissionDate)) %>%
  mutate(time_diff = difftime(submission_date, as.Date(collection_date), "weeks")) %>%
  separate(geo_loc_name, c("country"), sep = "\\:", remove = F)

morsels <- foreach(d = c("early", "alpha", "delta", "ba1", "ba5", "xbb", "pirola")) %do% {
  temp_df <- fread(str_glue("results/allele_frequency_out/codon_out/{d}.raw.filt.csv"))
}

codon_df <- bind_rows(morsels) %>%
  filter(id %in% names(fna)) %>%
  filter(freq < 0.5) %>%
  left_join(const_df %>% select(source_variant, mutation_name))

lineage_counts <- codon_df %>%
  group_by(id) %>%
  summarise(n_lineages = n_distinct(source_variant, na.rm = T),
            n_muts = n_distinct(mutation_name),
            n_var = sum(!is.na(source_variant)),
            sd_freq = sd(freq)) %>%
  left_join(meta %>% select(id = biosample, time_diff, country)) 

lineage_counts %>%
  left_join(meta %>% select(id = biosample, time_diff)) %>%
  ggplot(aes(x = time_diff, y = n_muts)) +
  geom_bin2d()

test <- glm(n_var ~ country + time_diff,
    data = lineage_counts,
    family = "poisson")

summary(test)
anova(test)
cor.test(lineage_counts$n_lineages, as.numeric(lineage_counts$time_diff), method = "spearman")
cor.test(lineage_counts$n_muts, as.numeric(lineage_counts$time_diff), method = "spearman")
cor.test(lineage_counts$n_var, as.numeric(lineage_counts$time_diff), method = "spearman")
cor.test(lineage_counts$sd_freq, as.numeric(lineage_counts$time_diff), method = "spearman")

median_time <- meta %>%
  filter(country != "") %>%
  group_by(country) %>% 
  summarise(median_lag = median(time_diff, na.rm = T)) %>%
  arrange(desc(median_lag))

country_plot_df <- meta %>%
  filter(country != "") %>%
  mutate(country = factor(country, unique(median_time$country))) 

pal <- distinctColorPalette(n_distinct(country_plot_df$country))

country_plot_df %>%
  ggplot(aes(x = country, y = time_diff, fill = country)) +
  geom_boxplot() +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Country", y = "Diff. between submission and collection dates")

meta %>%
  filter(time_diff > 1200) %>% View()
