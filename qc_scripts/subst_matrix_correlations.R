rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

data("BLOSUM45")
data("BLOSUM62")
data("BLOSUM100")
data("PAM30")
data("PAM250")
# 
# wag <- as.data.frame(fread("data/external_datasets/wag.dat", fill = T))
# rownames(wag) <- colnames(wag)
# # wag[upper.tri(wag)] <- t(wag)[upper.tri(wag)]

aa_list <- c("G", "A", "V", "L", "I",
             "T", "S", "M", "C", "P",
             "F", "Y", "W", "H", "K",
             "R", "D", "E", "N", "Q")

morsels <- foreach(mat_name = c("BLOSUM62", "BLOSUM45", "BLOSUM100", "PAM30", "PAM250")) %do% {
  mat <- get(mat_name)
  col_names <- colnames(mat)
  mat[upper.tri(mat, diag = T)] <- NA
  rownames(mat) <- col_names
  
  parsed <- as.data.frame(mat) %>%
    rownames_to_column("ref_AA") %>%
    pivot_longer(!ref_AA, names_to = "var_AA", values_to = "score") %>%
    filter(!is.na(score)) %>%
    mutate(matrix = mat_name)
}

merged <- bind_rows(morsels) %>%
  filter(ref_AA != var_AA) %>%
  filter(ref_AA %in% aa_list) %>%
  filter(var_AA %in% aa_list) %>% 
  pivot_wider(id_cols = c(ref_AA, var_AA), names_from = matrix, values_from = score) %>%
  select(-ref_AA, -var_AA)

as.data.frame(cor(merged)) %>%
  rownames_to_column("score1") %>%
  pivot_longer(!score1, names_to = "score2",  values_to = "score") %>%
  ggplot(aes(x = score1, y = score2, fill = score)) +
  geom_tile() +
  geom_text(aes(label = signif(score, 2))) +
  scale_fill_viridis_c() +
  labs(x = "Matrix", y = "Matrix", fill = "Pearson's R")
  
