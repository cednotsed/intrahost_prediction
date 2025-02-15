rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(randomcoloR)

fna <- readDNAStringSet("data/alignments/reassembled.masked.filt.aln")
df <- fread("data/metadata/all_sra_metadata.csv") %>%
  filter(biosample %in% names(fna))

plot_df <- df %>%
  separate(geo_loc_name, c("country"), "\\:") %>%
  mutate(country = ifelse(country == "", "Unknown", country))

pal <- distinctColorPalette(n_distinct(plot_df$country))

df %>%
  separate(geo_loc_name, c("country"), "\\:") %>%
  mutate(country = ifelse(country == "", "Unknown", country)) %>%
  group_by(dataset, country) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = dataset, y = n, fill = country)) +
  geom_col(color = "black") +
  scale_fill_manual(values = pal)

ggsave("results/qc_out/geographical_distribution.pdf", dpi = 600, width = 12, height = 8)
