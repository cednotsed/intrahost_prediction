rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(paletteer)

fna <- readDNAStringSet("data/alignments/reassembled.gap_filtered.pango_filtered.aln")
meta <- fread("data/metadata/all_sra_metadata.csv") %>%
  filter(biosample %in% names(fna))

library_counts <- meta %>%
  group_by(alias) %>%
  summarise(n_biosamples = n_distinct(biosample))

file_dir <- "results/mutation_stats/"
file_list <- list.files(file_dir, full.names = T)
file_list <- file_list[!grepl("_to_|spike|tonkin|dprime|max", file_list)]
file_list

morsels <- foreach(file_name = file_list) %do% {
  d <- gsub(file_dir, "", file_name)
  d <- gsub(".stats.csv", "", d)
  fread(file_name) %>%
    mutate(alias = d)
}

merged <- bind_rows(morsels) %>%
  as_tibble()

preds <- c("n", "max_freq", "median_freq", 
           "blosum62_score", "abs_charge", "abs_hydropathy",
           "abs_mw", "delta_charge", "delta_mw", 
           "delta_hydropathy", "delta_bind", "delta_expr", "mean_escape")

morsels <- foreach(d = unique(merged$alias)) %do% {
  foreach(pred = preds) %do% {
    temp <- merged %>%
      filter(alias == d)
    
    res <- cor.test(deframe(temp[, pred]), temp$global_n, method = "spearman")
    
    tibble(alias = d, predictor = pred, 
           rho = res$estimate, pval = res$p.value)
  }
}

plot_df <- bind_rows(morsels) %>% 
  mutate(adj_p = p.adjust(pval, method = "BH")) %>%
  mutate(not_significant = ifelse(adj_p > 0.05, "x", NA)) %>%
  mutate(alias = factor(alias, c("early", "alpha", "delta", 
                                 "ba1", "ba5", "xbb", 
                                 "pirola")))


order_df <- plot_df %>%
  group_by(predictor) %>%
  summarise(mean_rho = mean(rho)) %>%
  arrange(desc(mean_rho))

pal <- c(early = "#5A6F80FF", alpha = "#0E84B4FF", delta = "#B50A2AFF", 
         ba1 = "#E9D097FF", ba5 = "#278B9AFF", xbb = "#AE93BEFF", pirola = "grey70")

plot_df %>%
  mutate(predictor = factor(predictor, unique(order_df$predictor))) %>%
  ggplot(aes(x = predictor, y = rho, fill = alias)) +
  geom_col(position = "dodge",
           color = "black") +
  geom_text(aes(label = not_significant),
            position = position_dodge(width=0.9)) +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        legend.title = element_text(face = "bold")) +
  labs(x = "Predictor", y = "Spearman's rho", fill = "Dataset")
# ylim(-0.30, 0.8)

ggsave("results/prediction_out/predictor_spearman_correlation.pdf", dpi = 600, width = 5, height = 3)

