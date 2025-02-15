rm(list = ls())
setwd("C:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(randomcoloR)
require(ggpubr)

fna <- readDNAStringSet("data/alignments/tonkin.gap_filtered.aln")
meta <- fread("data/external_datasets/tonkin/tonkin.metadata.csv")
snp_df <- fread("results/pipeline_out.tonkin/snp_counts.csv")

pango_df <- fread("data/metadata/tonkin.gap_filtered.pango_lineage.csv") %>%
  dplyr::rename(sample_name = taxon) %>%
  left_join(snp_df)

pango_df %>%
  filter(conflict > 0) %>%
  group_by(lineage) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

plot_df <- pango_df %>%
  mutate(scorpio_call = ifelse(scorpio_call == "", "No designation", scorpio_call)) %>%
  group_by(scorpio_call) %>%
  mutate(outlier = is_outlier(n_SNPs))

length(fna) == nrow(plot_df)

pal <- distinctColorPalette(n_distinct(plot_df$scorpio_call))

before <- plot_df %>%
  ggplot() +
  geom_boxplot(aes(x = scorpio_call, y = n_SNPs, fill = scorpio_call)) +
  geom_text(aes(x = 2, y = 125, label = str_glue("n={n}")),
            data = before_counts) +
  facet_grid(cols = vars(dataset),
             scales = "free", space = "free") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  labs(x = "Scorpio lineage designation", y = "SNPs to Wuhan-Hu-1")

after <- plot_df %>%
  filter(!outlier) %>%
  ggplot() +
  geom_boxplot(aes(x = scorpio_call, y = n_SNPs, fill = scorpio_call)) +
  geom_text(aes(x = 2, y = 125, label = str_glue("n={n}")),
            data = after_counts) +
  facet_grid(cols = vars(dataset),
             scales = "free", space = "free") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  labs(x = "Scorpio lineage designation", y = "SNPs to Wuhan-Hu-1")


ggarrange(before, after, nrow = 2, align = "hv")
ggsave("results/qc_out/pango_filter.pdf", dpi = 600, height = 10, width = 12)

# Filter
plot_filt <- plot_df %>%
  filter(!outlier)

fna_filt <- fna[plot_filt$biosample]
writeXStringSet(fna_filt, "data/alignments/reassembled.gap_filtered.pango_filtered.aln")
