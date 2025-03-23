rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(paletteer)

file_dir <- "results/ML_out.020225/shap_out/"
file_list <- list.files(file_dir, full.names = T)
file_list <- file_list[!grepl("partition|spike|max|_to_|dprime\\.|linkage|prior", file_list)]

mut_meta <- fread("results/allele_frequency_out/observed_out/all_time.aggregate.csv") %>%
  distinct(mutation_name, protein_name, codon_number)

morsels <- foreach(file_name = file_list) %do% {
  id <- gsub(file_dir, "", file_name)
  id <- gsub(".shap.csv", "", id)
  
  fread(file_name) %>%
    mutate(alias = id)
}

merged <- bind_rows(morsels) 

shap_df <- merged %>%
  select(contains("shap"), alias, mutation_name) %>%
  pivot_longer(!c(alias, mutation_name), names_to = "predictor", values_to = "shap") %>%
  mutate(predictor = gsub(".shap", "", predictor))

value_df <- merged %>%
  select(!contains("shap"), alias) %>%
  pivot_longer(!c(alias, mutation_name), names_to = "predictor", values_to = "value")

parsed <- shap_df %>%
  bind_cols(value_df %>% select(value)) %>%
  ungroup() %>%
  left_join(mut_meta) %>%
  mutate(is_rbd = protein_name == "Spike" & codon_number %in% seq(331, 531)) %>%
  mutate(relevant = predictor %in% c("delta_bind", "delta_expr", "mean_escape"))

median_all_df <- parsed %>%
  group_by(predictor) %>%
  summarise(median_abs_shap = median(abs(shap))) %>%
  arrange(desc(median_abs_shap))

median_rbd_df <- parsed %>%
  filter(is_rbd) %>%
  group_by(predictor) %>%
  summarise(median_abs_shap = median(abs(shap))) %>%
  arrange(desc(median_abs_shap))

median_others_df <- parsed %>%
  filter(!is_rbd) %>%
  group_by(predictor) %>%
  summarise(median_abs_shap = median(abs(shap))) %>%
  arrange(desc(median_abs_shap))


pal <- c("grey", "#88C0D0FF")


# RBD
rbd <- parsed %>%
  filter(is_rbd) %>%
  mutate(predictor = factor(predictor, median_rbd_df$predictor)) %>%
  ggplot(aes(x = predictor, y = abs(shap), fill = relevant)) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Predictor", y = "Abs. SHAP values", fill = "DMS predictors")

others <- parsed %>%
  filter(!is_rbd) %>%
  mutate(predictor = factor(predictor, median_others_df$predictor)) %>%
  ggplot(aes(x = predictor, y = abs(shap), fill = relevant)) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Predictor", y = "Abs. SHAP values", fill = "DMS predictors")


ggarrange(rbd, others, nrow = 2, common.legend = T) 

ggsave("results/prediction_out/rbd_vs_non_rbd.shap.pdf", dpi = 600, height = 5, width = 8)
