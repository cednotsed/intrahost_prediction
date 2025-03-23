rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

possible_variants <- c("G", "A", "V", "L", "I",
                       "T", "S", "M", "C", "P",
                       "F", "Y", "W", "H", "K",
                       "R", "D", "E", "N", "Q")

meta <- fread("data/external_datasets/tonkin/tonkin.metadata.csv")

monthly_df <- fread("results/allele_frequency_out/observed_out/all_monthly_frequencies.csv") %>%
  mutate(collection_date = as.Date(str_glue("{collection_month}-01")))

timeframe <- meta %>%
    summarise(start = min(collection_month),
              end = max(collection_month)) %>%
    mutate(start = start %m-% months(1),
           end = end %m+% months(1)) %>%
    mutate(end_month = as.numeric(format(end, "%m")),
           end_year = as.numeric(format(end, "%Y")))
  
intra_stats <- fread(str_glue("results/mutation_stats/tonkin.stats.csv"))

monthly_filt <- monthly_df %>%
  filter(prop > 0.1) 

intra_filt <- intra_stats %>%
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
  mutate(month_diff = (12 * as.numeric(year) + as.numeric(month)) - (12 * timeframe$end_year + timeframe$end_month))

plt1 <- first_df %>% 
  filter(max_freq >= 0.5) %>%
  mutate(n_cat = case_when(n == 1 ~ "1",
                           n > 1 & n <= 10 ~ "2-10",
                           n > 10 ~ ">10")) %>%
  mutate(n_cat = factor(n_cat, c("1", "2-10", ">10"))) %>%
  ggplot(aes(x = month_diff, y = mutation_name, size = n_cat, color = max_freq)) + 
  geom_point() +
  labs(x = "Month until SAV reaches >10% freq.", y = "Intrahost SAV") +
  theme_bw() +
  scale_color_viridis_c()

plt2 <- first_df %>% 
  filter(max_freq < 0.5) %>%
  mutate(n_cat = case_when(n == 1 ~ "1",
                           n > 1 & n <= 10 ~ "2-10",
                           n > 10 ~ ">10")) %>%
  mutate(n_cat = factor(n_cat, c("1", "2-10", ">10"))) %>%
  ggplot(aes(x = month_diff, y = mutation_name, size = n_cat, color = max_freq)) + 
  geom_point() +
  labs(x = "Month until SAV reaches >10% freq.", y = "Intrahost SAV") +
  theme_bw() +
  scale_color_viridis_c()

ggarrange(plt1, plt2, nrow = 2,
          heights = c(2, 5))

ggsave("results/prediction_out/average_time_to_success.tonkin.pdf", dpi = 600, width = 8, height = 10)


first_df %>%
  filter(month_diff > 0) %>%
  filter(n > 1) %>%
  nrow()
