rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggrepel)
require(ggpubr)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")
meta <- fread("data/metadata/all_sra_metadata.csv") %>%
  distinct(alias)

morsels <- foreach(dataset = c("alpha", "delta", "ba1", "ba5", "xbb", "pirola")) %do% {
  # Intrahost
  intra_df <- fread(str_glue("results/linkage_out/intrahost_linkage.all/{dataset}.n5.with_zeroes.csv")) %>%
    mutate(mut1 = ifelse(mutation1 > mutation2, mutation1, mutation2),
           mut2 = ifelse(mutation1 > mutation2, mutation2, mutation1)) %>%
    dplyr::rename(intra_corr = corr)
  
  # Remove same codon mutations
  intra_hookup <- tibble(mutation_name = unique(c(intra_df$mut1, intra_df$mut2))) %>%
    separate(mutation_name, c("protein_name", "mut"), "_", remove = F) %>%
    mutate(codon_number = parse_number(mut)) %>%
    mutate(codon_name = str_glue("{protein_name}_{codon_number}")) %>%
    select(protein_name, mutation_name, codon_name)
  
  intra_filt <- intra_df %>% 
    left_join(intra_hookup %>% select(mut1 = mutation_name, codon1 = codon_name, protein_name1 = protein_name)) %>%
    left_join(intra_hookup %>% select(mut2 = mutation_name, codon2 = codon_name, protein_name2 = protein_name)) %>%
    filter(codon1 != codon2) %>%
    mutate(within_protein = protein_name1 == protein_name2) %>%
    mutate(is_linked = intra_corr > 0.9)
  
  intra_mat <- intra_filt %>%
    group_by(is_linked, within_protein) %>%
    summarise(n = n()) %>%
    select(is_linked, within_protein, n) %>%
    pivot_wider(id_cols = is_linked, names_from = within_protein, values_from = n) %>%
    column_to_rownames("is_linked")
  
  intra_test <- fisher.test(intra_mat)
  
  # Consensus
  consensus_df <- fread(str_glue("results/linkage_out/observed_linkage/observed_Dprime.before_{dataset}.gt1000.csv"))
  
  consensus_hookup <- tibble(mutation_name = unique(c(consensus_df$mut1, consensus_df$mut2))) %>%
    separate(mutation_name, c("protein_name", "mut"), "_", remove = F) %>%
    mutate(codon_number = parse_number(mut)) %>%
    mutate(codon_name = str_glue("{protein_name}_{codon_number}")) %>%
    select(protein_name, mutation_name, codon_name)
  
  consensus_filt <- consensus_df %>% 
    left_join(consensus_hookup %>% select(mut1 = mutation_name, codon1 = codon_name, protein_name1 = protein_name)) %>%
    left_join(consensus_hookup %>% select(mut2 = mutation_name, codon2 = codon_name, protein_name2 = protein_name)) %>%
    filter(codon1 != codon2) %>%
    mutate(within_protein = protein_name1 == protein_name2) %>%
    mutate(is_linked = Dprime > 0.9)
  
  consensus_mat <- consensus_filt %>%
    group_by(is_linked, within_protein) %>%
    summarise(n = n()) %>%
    select(is_linked, within_protein, n) %>%
    pivot_wider(id_cols = is_linked, names_from = within_protein, values_from = n) %>%
    column_to_rownames("is_linked")
  
  consensus_test <- fisher.test(consensus_mat)
  
  merged <- consensus_filt %>%
    inner_join(intra_filt)
  
  corr1 <- cor.test(merged$corr_no_double_zeroes, merged$Dprime, method = "pearson")
  r1 <- signif(corr1$estimate, 2)
  pval1 <- signif(corr1$p.value, 2)
  n_consensus <- nrow(consensus_filt)
  n_intra <- nrow(intra_filt)
  n_pairs <- nrow(merged)
  n_intra_linked <- sum(intra_filt$intra_corr > 0.9)
  n_consensus_linked <- sum(consensus_filt$Dprime > 0.9)
  
  plt <- merged %>% 
    ggplot(aes(x = intra_corr, y = Dprime)) +
    geom_bin2d() +
    geom_smooth(method = "lm") +
    theme_bw() +
    scale_fill_viridis_c(trans = "log10") + 
    theme(text = element_text(family = "sans"),
          axis.title = element_text(face = "bold"),
          legend.position = "none") +
    labs(x = "Intrahost linkage (r)", y = "Consensus linkage (D')",
         title = str_glue("dataset={dataset}, r={r1}, n={n_pairs}")) +
    ylim(-1, 1) +
    xlim(-1, 1)
  
  # return(tibble(alias = dataset, n_pairs = n_pairs, n_intra = n_intra,
  #               n_intra_linked = n_intra_linked))
  return(plt)
}

ggarrange(plotlist = morsels)

ggsave("results/linkage_out/all_intra_versus_cons_heatmaps.with_labels.pdf", dpi = 600, width = 12, height = 8)
morsels[[2]] 

ggsave("results/linkage_out/delta_intra_versus_cons_heatmap.pdf", dpi = 600, width = 5, height = 3)
