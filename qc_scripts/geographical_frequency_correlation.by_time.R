rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)

agg_df <- fread("results/allele_frequency_out/observed_out/all_time.aggregate.csv")

# Get monthly chunks
month_list <- deframe(fread("results/allele_frequency_out/observed_out/all_monthly_frequencies.csv") %>%
  distinct(collection_month) %>%
    mutate(collection_month = as.Date(collection_month)))

chunks <- split(month_list, ceiling(seq_along(month_list)/ 8))

chunk_meta <- foreach(i = seq(length(chunks)), .combine = "bind_rows") %do% {
  chunk <- chunks[[i]]
  
  tibble(collection_month = chunk, chunk_index = i) %>%
    mutate(start = min(collection_month),
           end = max(collection_month)) %>%
    mutate(start = format(start, "%b-%y"),
           end = format(end, "%b-%y"))
}

region_list <- c("asia", "oceania", "europe", "north_america", "south_america", "africa")

pairs <- combn(region_list, m = 2)

time_morsels <- foreach(idx = seq(ncol(pairs))) %do% {
  pair <- pairs[, idx]
  region1 <- pair[1]
  region2 <- pair[2]
  
  df1 <- fread(str_glue("results/allele_frequency_out/observed_out/geographical_splits/{region1}.all_monthly_frequencies.csv")) %>%
    filter(mutation_name %in% agg_df$mutation_name) %>%
    mutate(collection_month = as.Date(collection_month)) %>%
    left_join(chunk_meta)
  
  df2 <- fread(str_glue("results/allele_frequency_out/observed_out/geographical_splits/{region2}.all_monthly_frequencies.csv")) %>%
    filter(mutation_name %in% agg_df$mutation_name) %>%
    mutate(collection_month = as.Date(collection_month)) %>%
    left_join(chunk_meta)
  
  # calculate r by chunk
  df1_agg <- df1 %>%
    group_by(mutation_name, chunk_index, start, end) %>%
    summarise(prop = sum(n_present) / sum(n_total))
  
  df2_agg <- df2 %>%
    group_by(mutation_name, chunk_index, start, end) %>%
    summarise(prop = sum(n_present) / sum(n_total))
  
  name1 <- Hmisc::capitalize(region1)
  name2 <- Hmisc::capitalize(region2)
  
  plot_df <- df1_agg %>%
    dplyr::rename(prop1 = prop) %>%
    full_join(df2_agg %>% dplyr::rename(prop2 = prop)) %>%
    mutate(prop1 = replace_na(prop1, 0),
           prop2 = replace_na(prop2, 0)) %>%
    group_by(chunk_index) %>%
    summarise(r = cor(prop1, prop2)) %>%
    left_join(chunk_meta %>%
                distinct(chunk_index, start, end)) %>%
    mutate(chunk_name = str_glue("{start}->{end}")) %>%
    mutate(comparison = str_glue("{name1}-{name2}"))
  
  return(plot_df)
  # 
  # plot_df %>%
  #   mutate(chunk_name = factor(chunk_name, unique(plot_df$chunk_name))) %>%
  #   ggplot(aes(x = chunk_name, y = r)) +
  #   geom_col(color = "black") +
  #   theme_bw() +
  #   ylim(0, 1) +
  #   labs(x = str_glue("{name1}"), y = str_glue("{name2}"))
}

combined <- bind_rows(time_morsels)

combined %>%
  mutate(chunk_name = factor(chunk_name, unique(combined$chunk_name))) %>%
  ggplot(aes(x = chunk_name, y = r, fill = comparison)) + 
  geom_col(color = "black") +
  facet_wrap(comparison ~ ., nrow = 3, ncol = 5) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Time period", y = "Pearson's r")

ggsave("results/qc_out/geographical_correlation.by_time.pdf", dpi = 300, width = 10, height = 5)


# merged <- df1 %>%
#   select(mutation_name, collection_month, region1 = prop) %>%
#   full_join(df2 %>% 
#               select(mutation_name, collection_month, region2 = prop)) %>%
#   ungroup() %>%
#   complete(collection_month, mutation_name) %>%
#   mutate(region1 = replace_na(region1, 0),
#          region2 = replace_na(region2, 0)) %>%
#   mutate(diff = region1 - region2)
# 
# merged %>%
#   ggplot(aes(x = collection_month, y = diff)) +
#   geom_bin2d(bins = 100)
# 
# # # Check all months present
# # merged %>%
# #   group_by(mutation_name) %>%
# #   summarise(n = n()) %>%
# #   distinct(n)
# 
# 
# 
