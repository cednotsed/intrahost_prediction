rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(paletteer)

file_dir <- "results/ML_out/shap_out/"
file_list <- list.files(file_dir, full.names = T)
file_list <- file_list[!grepl("_to_|dprime\\.|linkage|prior", file_list)]

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
  filter(grepl("partition.dprime_only", alias)) %>%
  ungroup()

parsed %>% distinct(alias)

parsed %>%
  group_by(predictor) %>%
  summarise(median_abs_shap = median(abs(shap))) %>%
  arrange(desc(median_abs_shap))

# Add mutation annotations
dat_morsels <- foreach(d = c("alpha", "delta", "ba1", "ba5", "xbb", "pirola")) %do% {
  stat_df <- fread(str_glue("results/prediction_out/mutation_stats/{d}.dprime_stats.csv")) %>%
    filter(max_dprime != -1)
  
  plot_temp <- parsed %>%
    filter(alias == str_glue("{d}.partition.dprime_only")) %>%
    right_join(stat_df %>%
                select(mutation_name, future_high))

  return(plot_temp)  
}

intra_corr <- c("n_linked_binary", "n_linked", "max_corr")
consensus_corr <- c("n_dprime_linked_binary", "n_dprime_linked", "max_dprime")
intrahost <- c("n", "max_freq", "median_freq")
physio <- c("blosum62_score", "delta_mw", "delta_charge",
            "delta_hydropathy", "abs_mw", "abs_charge",
            "abs_hydropathy")
dms <- c("delta_bind", "delta_expr", "mean_escape")

plot_df <- bind_rows(dat_morsels) %>% 
  filter(predictor %in% c(intra_corr, consensus_corr)) %>% 
  group_by(mutation_name, alias, future_high) %>%
  summarise(sum_shap = sum(shap)) %>%
  mutate(alias = gsub(".partition.dprime_only", "", alias)) %>%
  mutate(alias = factor(alias, c("alpha", "delta", "ba1", "ba5", "xbb", "pirola")))

count_df <- plot_df %>%
  group_by(alias, future_high) %>%
  summarise(n = n_distinct(mutation_name))

plot_df %>%
  ggplot(aes(x = alias, y = sum_shap, fill = future_high)) +
  geom_boxplot() +
  geom_text(aes(y = 3, label = str_glue("{n}")),
            data = count_df,
            position = position_dodge(width = .9),
            size = 2) +
  geom_pwc(label.size = 2) +
  theme_bw() +
  scale_fill_manual(values = c("indianred", "royalblue")) +
  theme(text = element_text(family = "sans"),
        axis.title = element_text(face = "bold"),
        legend.position = "none") +
  labs(x = "Dataset", y = "Impact of linkage predictors (SHAP)")

ggsave("results/prediction_out/linkage_impact.all.pdf", dpi = 600, height = 3, width = 4.5)
fwrite(plot_df, "results/prediction_out/linkage_impact.csv")
