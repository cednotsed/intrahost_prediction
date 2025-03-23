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

morsels <- foreach(file_name = file_list) %do% {
  id <- gsub(file_dir, "", file_name)
  id <- gsub(".shap.csv", "", id)
  
  fread(file_name) %>%
    mutate(alias = id)
}

merged <- bind_rows(morsels) %>%
  filter(!grepl("_to_|spike|partition|max", alias)) %>%
  select(contains("shap"), alias) %>%
  pivot_longer(!alias, names_to = "predictor", values_to = "shap") %>%
  filter(!is.na(shap))

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

ggsave("results/prediction_out/predictor_importance.pdf", dpi = 600, width = 4.5, height = 3)
  
  
