rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(doParallel)

meta <- fread("data/metadata/sra_metadata/filtered_sra_accessions.transition.csv") %>%
  select(id = biosample, collection_date) %>%
  mutate(collection_month = format(collection_date, "%Y-%m"))

monthly_df <- fread("results/allele_frequency_out/all_monthly_frequencies.csv") %>%
  mutate(collection_date = as.Date(paste0(collection_month, "-01")))

monthly_agg <- monthly_df %>%
  group_by(mutation_name, collection_month, collection_date) %>%
  summarise(monthly_n = sum(n_present),
            monthly_total = unique(n_total)) %>%
  mutate(monthly_prop = monthly_n / monthly_total)

current_agg <- monthly_agg %>%
  filter(collection_month == "2021-06") %>%
  filter(monthly_prop > 0.1)
  
parsed_agg <- monthly_agg %>%
  filter(collection_date > as.Date("2022-01-01")) %>%
  group_by(mutation_name) %>%
  summarise(global_n = sum(monthly_n),
            max_prop = max(monthly_prop))

parsed_filt <- parsed_agg %>%
  filter(max_prop > 0.1) %>%
  filter(!(mutation_name %in% current_agg$mutation_name))

freq_df <- fread("results/allele_frequency_out/all_sites.transition.raw.csv")
agg_df <- fread("results/allele_frequency_out/all_sites.transition.csv")

detect_df <- freq_df %>%
  left_join(meta) %>%
  filter(collection_month == "2021-06")
filter(mutation_name %in% parsed_filt$mutation_name) %>% 
  distinct(mutation_name)

plot_df <- agg_df %>%
  mutate(is_fixed = mutation_name %in% parsed_filt$mutation_name)
  
plot_df %>%
  ggplot(aes(x = is_fixed, y = n)) +
  geom_boxplot()

logreg <- glm(is_fixed ~ n + median_freq + max_freq,
              )


# detect_df <- agg_df %>%
#   # filter(max_freq < 0.9) %>% 
#   filter(collection_month == "2021-06") %>%
#   filter(mutation_name %in% parsed_filt$mutation_name) %>% 
#   left_join(parsed_agg)

monthly_df %>%
  filter(mutation_name %in% detect_df$mutation_name) %>%
  ggplot(aes(x = collection_date, y = prop, color = mutation_name)) +
  annotate("rect", 
           xmin = as.Date("2021-06-01"), 
           xmax = as.Date("2021-12-01"), 
           ymin = 1, 
           ymax = 1.2, 
           fill = "#9E8356FF",
           alpha = 0.5) +
  annotate("rect", 
           xmin = as.Date("2021-12-01"), 
           xmax = as.Date("2023-02-01"), 
           ymin = 1, 
           ymax = 1.2, 
           fill = "#278B9AFF",
           alpha = 0.5) +
  geom_point() + 
  geom_line() +
  geom_vline(xintercept = as.Date("2021-06-01")) +
  geom_hline(yintercept = 0.1, lty = "dashed") +
  scale_x_date(date_breaks = "2 months", 
               date_labels = "%b-%y",
               limits = c(as.Date("2019-12-01"), as.Date("2024-07-01"))) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(family = "sans"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  labs(x = "Collection month", y = "Prop. genomes", color = "Mutation")

parsed_filt %>%
  filter(mutation_name %in% detect_df$mutation_name) %>%
  left
