rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)

agg_df <- fread("results/allele_frequency_out/observed_out/all_time.aggregate.csv")

region_list <- c("asia", "oceania", "europe", "north_america", "south_america", "africa")

pairs <- combn(region_list, m = 2)

overall_morsels <- foreach(idx = seq(ncol(pairs))) %do% {
  pair <- pairs[, idx]
  region1 <- pair[1]
  region2 <- pair[2]
  
  df1 <- fread(str_glue("results/allele_frequency_out/observed_out/geographical_splits/{region1}.all_monthly_frequencies.csv")) %>%
    filter(mutation_name %in% agg_df$mutation_name)

  df2 <- fread(str_glue("results/allele_frequency_out/observed_out/geographical_splits/{region2}.all_monthly_frequencies.csv")) %>%
    filter(mutation_name %in% agg_df$mutation_name)

  # All time
  total_df1 <- df1 %>%
    group_by(mutation_name) %>%
    summarise(prop = sum(n_present) / sum(n_total))
  
  total_df2 <- df2 %>%
    group_by(mutation_name) %>%
    summarise(prop = sum(n_present) / sum(n_total))
  
  prop_agg <- total_df1 %>%
    dplyr::rename(prop1 = prop) %>%
    full_join(total_df2 %>% dplyr::rename(prop2 = prop)) %>%
    mutate(prop1 = replace_na(prop1, 0),
           prop2 = replace_na(prop2, 0)) 
  
  r <- signif(cor(prop_agg$prop1, prop_agg$prop2), 2)
  
  name1 <- Hmisc::capitalize(region1)
  name2 <- Hmisc::capitalize(region2)
  
  plt <- prop_agg %>%
    ggplot(aes(x = prop1, y = prop2)) +
    geom_point(color = "orange", alpha = 0.5) +
    geom_smooth(method = "lm") +
    annotate("text", x = 0, y = 1, 
             label = str_glue("r={r}"),
             hjust = 0) +
    theme_bw() +
    labs(x = str_glue("{name1}"), y = str_glue("{name2}"))
  
  return(plt)
}

combined <- ggarrange(plotlist = overall_morsels, nrow = 3, ncol = 5)

ggsave("results/qc_out/geographical_correlation.overall.png", dpi = 300, width = 10, height = 5)


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
