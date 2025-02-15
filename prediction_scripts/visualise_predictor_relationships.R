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

morsels <- foreach(file_name = file_list) %do% {
  id <- gsub(file_dir, "", file_name)
  id <- gsub(".shap.csv", "", id)
  
  fread(file_name) %>%
    mutate(alias = id)
}

set.seed(66)

merged <- bind_rows(morsels) %>%
  filter(!grepl("_to_|linkage", alias)) %>%
  filter(grepl("spike", alias)) %>%
  # group_by(alias) %>%
  # sample_n(50, replace = F) %>%
  ungroup()

shap_df <- merged %>%
  select(contains("shap"), alias, mutation_name) %>%
  pivot_longer(!c(alias, mutation_name), names_to = "predictor", values_to = "shap") %>%
  mutate(predictor = gsub(".shap", "", predictor))

value_df <- merged %>%
  select(!contains("shap"), alias) %>%
  pivot_longer(!c(alias, mutation_name), names_to = "predictor", values_to = "value")

plot_df <- shap_df %>%
  bind_cols(value_df %>% select(value)) %>%
  filter(!is.na(shap))

plot_df %>%
  summarise(n_distinct(predictor),
            n_distinct(alias))

preds <- deframe(plot_df %>% distinct(predictor))

plot_list <- list()
i <- 0

foreach(pred = preds) %do% {
  temp1 <- plot_df %>%
    filter(predictor == pred)
  
  if(pred %in% c("delta_bind", "delta_expr")) {
    temp1 <- temp1 %>%
      filter(value != -100)
  } else if (pred == "mean_escape") {
    temp1 <- temp1 %>%
      filter(value != -1)
  }
  
  crumbs <- foreach(d = c("early.spike", "alpha.spike", "delta.spike", "ba1.spike", "ba5.spike", "xbb.spike", "pirola.spike")) %do% {
    i <- i + 1
    temp2 <- temp1 %>%
      filter(alias == d)
    
    median_value <- deframe(temp2 %>%
                              ungroup() %>%
                              summarise((min(value) + max(value)) / 2))
    
    plt <- temp2 %>%
      # sample_n(500, replace = F) %>%
      ggplot(aes(x = value, y = shap)) +
      geom_point(color = "steelblue") +
      geom_hline(yintercept = 0, 
                 lty = "dashed", 
                 color = "black",
                 size = 1.5) +
      geom_smooth(aes(x = value, y = shap),
                  method = "loess", color = "indianred4",
                  data = temp2,
                  size = 2) +
      theme_bw() +
      # theme(axis.text = element_blank(),
      #       axis.title = element_blank(),
      #       axis.ticks = element_blank(),
      #       panel.grid = element_blank(),
      #       legend.position = "none")
      theme(text = element_text(family = "sans"),
            axis.title = element_text(face = "bold"),
            legend.position = "none") +
      labs(x = pred, y = "SHAP value", title = d)
    
    plot_list <- c(plot_list, list(plt))
    return(NULL)
  }
  
  return(NULL)
}

merged_plt <- ggarrange(plotlist = plot_list, align = "hv", ncol = 7, nrow = 13)
# ggsave("results/ML_out/all_predictor_relationships.png", plot = merged_plt, height = 20, width = 20)
ggsave("results/ML_out/all_predictor_relationships.with_labels.png", plot = merged_plt, height = 20, width = 20)


pred_list <- rev(c("max_freq", "delta_hydropathy", "blosum62_score", "delta_bind"))

plot_list2 <- foreach(pred = pred_list) %do% {
  plot_df %>%
    filter(alias == "delta.spike") %>%
    filter(predictor == pred) %>%
    filter(value != -100) %>%
    ggplot(aes(x = value, y = shap)) +
    geom_point(color = "darkseagreen4") +
    geom_hline(yintercept = 0, lty = "dashed", color = "grey") +
    geom_smooth(method = "loess", color = "blue") +
    theme_classic() +
    theme(text = element_text(family = "sans"),
          axis.title = element_text(face = "bold"),
          legend.position = "none") +
    labs(x = pred, y = "SHAP value")
}

ggarrange(plotlist = plot_list2, align = "hv", nrow = 1)
ggsave("results/ML_out/predictors_of_interest.pdf", height = 3, width = 12)
