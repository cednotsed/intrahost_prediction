rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)
require(see)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

possible_variants <- c("G", "A", "V", "L", "I",
                       "T", "S", "M", "C", "P",
                       "F", "Y", "W", "H", "K",
                       "R", "D", "E", "N", "Q")

meta <- fread("data/metadata/all_sra_metadata.csv")

d <- "delta"

monthly_df <- fread("results/allele_frequency_out/observed_out/all_monthly_frequencies.csv") %>%
  mutate(collection_date = as.Date(str_glue("{collection_month}-01")))

monthly_filt <- monthly_df %>%
  filter(prop > 0.1) 

time_morsels <- foreach(d = unique(meta$alias)) %do% {
  timeframes <- meta %>%
    group_by(dataset, alias) %>%
    summarise(start = min(collection_month),
              end = max(collection_month)) %>%
    mutate(start = start %m-% months(1),
           end = end %m+% months(1)) %>%
    mutate(end_month = as.numeric(format(end, "%m")),
           end_year = as.numeric(format(end, "%Y")))
  
  timeframe <- timeframes %>%
    filter(alias == d)
  
  mutation_stats <- fread(str_glue("results/mutation_stats/{d}.stats.csv")) %>%
    filter(future_high)
  
  intra_filt <- mutation_stats %>%
    filter(mutation_name %in% unique(monthly_filt$mutation_name))
  
  first_df <- monthly_filt %>%
    # filter(collection_date > timeframe$end) %>% 
    filter(mutation_name %in% intra_filt$mutation_name) %>%
    group_by(mutation_name) %>%
    arrange(collection_month) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    left_join(intra_filt %>% select(mutation_name, n, max_freq)) %>%
    separate(collection_month, c("year", "month"), "-") %>%
    mutate(month_diff = (12 * as.numeric(year) + as.numeric(month)) - (12 * timeframe$end_year + timeframe$end_month)) %>%
    mutate(alias = d)
  
  # bind_rows(morsels) %>%
  #   select(mutation_name, collection_month) %>%
  #   separate(collection_month, c("year", "month"), "-") %>%
  #   mutate(month_diff = (12 * as.numeric(year) + as.numeric(month)) - (12 * timeframe$end_year + timeframe$end_month)) %>%
  #   mutate(alias = d)
  
  return(first_df)
}

merged <- bind_rows(time_morsels)

merged %>%
  group_by(alias) %>%
  summarise(median_lag = median(month_diff)) %>%
  arrange(desc(median_lag))

pal <- c("#5A6F80FF", "#0E84B4FF", "#B50A2AFF","#E9D097FF", "#278B9AFF", "#AE93BEFF", "grey30")

plot_df <- merged %>%
  mutate(variant_annot = case_when(grepl("early", alias) ~ "Early",
                                   grepl("alpha", alias) ~ "Alpha",
                                   grepl("delta", alias) ~ "Delta",
                                   grepl("ba1", alias) ~ "BA.1",
                                   grepl("ba5", alias) ~ "BA.5",
                                   grepl("xbb", alias) ~ "XBB",
                                   grepl("pirola", alias) ~ "BA.2.86")) %>%
  mutate(variant_annot = factor(variant_annot, c("Early", "Alpha", "Delta", "BA.1", 
                                                 "BA.5", "XBB", "BA.2.86"))) %>%
  mutate(type = ifelse(max_freq >= 0.5, "Consensus SAVs", "Subconsensus SAVs"))

mut_counts <- plot_df %>%
  group_by(variant_annot, type) %>%
  summarise(n = n_distinct(mutation_name))

plot_df %>%  
  ggplot(aes(x = variant_annot, y = month_diff, fill = variant_annot)) +
  geom_point(color = "darkgrey", position = position_jitter(width = 0.08), 
             size = 0.5, 
             alpha = 0.8) +
  geom_boxplot(position = position_nudge(x = 0.3, y = 0), 
               width = 0.2, 
               outlier.shape = NA) +
  geom_text(aes(x = variant_annot, y = 0, label = str_glue("n={n}")),
            data = mut_counts) +
  facet_grid(rows = vars(type)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = pal) +
  labs(y = "Months to frequency peak", x = "Timeframe")

ggsave("results/prediction_out/average_time_to_success.all.pdf", dpi = 600, width = 7, height = 5)
