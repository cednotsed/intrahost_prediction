rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(paletteer)

dms_df <- fread("data/external_datasets/dms/parsed_dms_phenotypes.csv")

file_dir <- "results/ML_out.020225/shap_out/"
file_list <- list.files(file_dir, full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  id <- gsub(file_dir, "", file_name)
  id <- gsub(".shap.csv", "", id)
  
  fread(file_name) %>%
    mutate(alias = id)
}

merged <- bind_rows(morsels) %>%
  filter(alias %in% c("early", "alpha", "delta", "ba1", "ba5", "xbb", "pirola")) %>%
  select(mutation_name, contains("shap"), alias) %>%
  pivot_longer(!c(mutation_name, alias), names_to = "predictor", values_to = "shap") %>%
  filter(!is.na(shap)) %>%
  mutate(is_RBD = mutation_name %in% dms_df$mutation_name)

plot_df <- merged %>%
  filter(grepl("bind|escape|expr|max_freq", predictor)) %>%
  mutate(predictor = case_when(predictor == "delta_bind.shap" ~ "RBD binding",
                               predictor == "delta_expr.shap" ~ "RBD expression",
                               predictor == "mean_escape.shap" ~ "Abx. escape",
                               predictor == "max_freq.shap" ~ "Max. intrahost freq.")) 

order_df <- plot_df %>%
  filter(is_RBD) %>%
  group_by(predictor) %>%
  summarise(median = median(abs(shap))) %>%
  arrange(desc(median))

plot_df %>%
  mutate(predictor = factor(predictor, order_df$predictor)) %>%
  ggplot(aes(x = predictor, y = abs(shap), fill = is_RBD)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic() +
  geom_pwc() +
  labs(x = "Predictor", y = "Abs. SHAP value", fill = "Is RBD?")

median_df <- merged %>%
  group_by(predictor) %>%
  summarise(median_abs_shap = median(abs(shap))) %>%
  arrange(desc(median_abs_shap))

pal <- c(paletteer::paletteer_d("beyonce::X10")[2:5], paletteer::paletteer_d("nord::frost"),
         paletteer::paletteer_d("nord::lumina"))

merged %>%
  mutate(abs_shap = abs(shap)) %>%
  mutate(predictor = factor(predictor, median_df$predictor)) %>%
  ggplot(aes(x = predictor, y = abs_shap, fill = predictor)) +
  geom_boxplot(outliers = FALSE) +
  theme_bw() +
  scale_fill_manual(values = pal) +
  theme(text = element_text(family = "sans"),
        legend.position = "none",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Predictor", y = "Abs. SHAP values")

# ggsave("results/ML_out/predictor_importance.pdf", dpi = 600, width = 4.5, height = 3)

df <- fread("results/mutation_stats/alpha.stats.csv")

df %>% View()
