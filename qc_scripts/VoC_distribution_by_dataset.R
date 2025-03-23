rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(randomcoloR)
require(paletteer)

fna <- c(readDNAStringSet("data/alignments/reassembled.gap_filtered.pango_filtered.aln"))

pango_df <- fread("data/metadata/reassembled.gap_filtered.pango_lineage.csv") %>%
  select(biosample = taxon, lineage, conflict, scorpio_call, scorpio_support, scorpio_conflict)

meta <- fread("data/metadata/all_sra_metadata.csv") %>%
  mutate(collection_month = format(as.Date(collection_date, "%Y-%m-%d"), "%Y-%m"))

merged_meta <- pango_df %>%
  filter(biosample %in% names(fna)) %>%
  left_join(meta %>% select(biosample, bioproject, experiment, 
                            run, collection_date, collection_month,
                            geo_loc_name, dataset)) %>%
  mutate(collection_month = ifelse(is.na(collection_month), collection_date, collection_month)) %>%
  arrange(collection_month)

merged_meta %>% View()

plot_df <- merged_meta %>%
  mutate(voc = case_when(grepl("XBB", scorpio_call) ~ "XBB.*",
                         grepl("Alpha|B.1.1.7", scorpio_call) ~ "Alpha",
                         grepl("Delta", scorpio_call) ~ "Delta",
                         grepl("BA.1", scorpio_call) ~ "Omicron (BA.1)",
                         grepl("BA.5", scorpio_call) ~ "Omicron (BA.5)",
                         grepl("Beta|Epsilon|Iota|Lambda|Mu|Gamma|Eta|Zeta|B.1.1.318-like|XE|A.23|AZ", scorpio_call) ~ "Other lineages",
                         grepl("BA.2.|BA.4|Unassigned", scorpio_call) ~ "Omicron (others)",
                         scorpio_call == "" ~ "No designation")) %>%
  mutate(voc = ifelse(grepl("BA.2.86|JN.1", lineage), "BA.2.86/JN.1", voc)) %>%
  group_by(dataset, voc) %>%
  summarise(n = n_distinct(biosample))

total_counts <- merged_meta %>%
  group_by(dataset) %>%
  summarise(n_total = n())

pal <- distinctColorPalette(n_distinct(plot_df$voc))
pal <- c("#D56274", "#DBBD6D", "#D7D6C8",
         "#CB53DA", "#84B5D0", "#ABE45D",
         "#89E0B6", "#7D76D1", "#D9A4D2")

pal <- c("#009392FF", "#72AAA1FF", "#B1C7B3FF", 
         "#F1EAC8FF", "#E5B9ADFF", "#D98994FF", 
         "#5F93ACFF", "#778868FF", "#3D5941FF")

plot_df %>%
  mutate(dataset = factor(dataset, unique(merged_meta$dataset))) %>%
  ggplot(aes(x = dataset, y = n)) +
  geom_col(aes(fill = voc), color = "black") +
  geom_text(aes(x = dataset, y = n_total, label = str_glue("n={n_total}")),
            vjust = -0.1,
            data = total_counts) +
  scale_fill_manual(values = pal) +
  ylim(0, 1400) +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(family = "sans"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  labs(y = "No. biosamples", x = "Collection month", fill = "Lineage")

ggsave("results/qc_out/VoC_distribution_by_dataset.png", 
       dpi = 600,
       height = 4,
       width = 8)
