rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

possible_variants <- c("G", "A", "V", "L", "I",
                       "T", "S", "M", "C", "P",
                       "F", "Y", "W", "H", "K",
                       "R", "D", "E", "N", "Q")

meta <- fread("data/metadata/all_sra_metadata.csv")

# dataset <- "delta"

for(dataset in unique(meta$alias)) {
  intra_agg <- fread(str_glue("results/allele_frequency_out/codon_out/{dataset}.csv")) %>%
    filter(ref_AA %in% possible_variants) %>%
    filter(var_AA %in% possible_variants)
  
  # Prior frequency
  before_agg <- fread(str_glue("results/allele_frequency_out/observed_out/before_{dataset}.aggregate.csv")) %>%
    select(mutation_name, prior_n = global_n)
  
  # High freq mutations during dataset timeframe
  current_agg <- fread(str_glue("results/allele_frequency_out/observed_out/during_{dataset}.aggregate.csv"))
  
  # Future frequency
  future_agg <- fread(str_glue("results/allele_frequency_out/observed_out/after_{dataset}.aggregate.csv"))
  
  merged <- intra_agg %>%
    left_join(future_agg %>% dplyr::select(mutation_name, max_prop, global_n, global_prop)) %>%
    left_join(before_agg) %>%
    mutate(prior_n = replace_na(prior_n, 0)) %>%
    mutate(global_n = replace_na(global_n, 0)) %>%
    mutate(global_prop = replace_na(global_prop, 0)) %>%
    mutate(max_prop = replace_na(max_prop, 0)) %>%
    mutate(future_high = max_prop > 0.1) %>%
    mutate(future_fixed = max_prop > 0.9) %>% 
    filter(!(mutation_name %in% current_agg$mutation_name)) # Remove mutations currently at high frequency
  
  merged %>% 
    fwrite(str_glue("results/mutation_stats/{dataset}.stats.csv"),
           eol = "\n")
  
  table(merged$future_fixed)
  table(merged$future_high)
}
  # p1 <- merged %>%
  #   ggplot(aes(y = log10(global_n + 1), x = max_freq)) +
  #   geom_bin2d() +
  #   scale_fill_viridis_c(trans = "log10") +
  #   geom_smooth(method = "lm", color = "black", fill = "grey") +
  #   theme_bw() +
  #   theme(legend.position = "top",
  #         text = element_text(family = "sans"),
  #         axis.title = element_text(face = "bold"),
  #         legend.title = element_text(face = "bold")) +
  #   labs(x = "Max intrahost freq.", y = "Log10(no. GISAID sequences)",
  #        fill = "No. mutations") 
  # 
  # p2 <- merged %>%
  #   ggplot(aes(y = log10(global_n + 1), x = n)) +
  #   geom_bin2d() +
  #   scale_fill_viridis_c(trans = "log10") +
  #   geom_smooth(method = "lm", color = "black", fill = "grey") +
  #   theme_bw() +
  #   theme(legend.position = "top",
  #         text = element_text(family = "sans"),
  #         axis.title = element_text(face = "bold"),
  #         legend.title = element_text(face = "bold")) +
  #   labs(x = "No. libraries detected", y = "Log10(no. GISAID sequences)",
  #        fill = "No. mutations") 
  # 
  # p3 <- merged %>%
  #   ggplot(aes(x = future_high, y = log10(max_freq))) +
  #   geom_boxplot() +
  #   geom_pwc() +
  #   geom_boxplot(fill = "turquoise") +
  #   geom_pwc() +
  #   theme_bw() +
  #   theme(legend.position = "top",
  #         text = element_text(family = "sans"),
  #         axis.title = element_text(face = "bold"),
  #         legend.title = element_text(face = "bold")) +
  #   labs(x = "Max. future frequency >10%", 
  #        y = "Max. intrahost freq.")
  # 
  # p4 <- merged %>%
  #   ggplot(aes(x = future_high, y = log10(n))) +
  #   geom_boxplot(fill = "goldenrod") +
  #   geom_pwc() +
  #   theme_bw() +
  #   theme(legend.position = "top",
  #         text = element_text(family = "sans"),
  #         axis.title = element_text(face = "bold"),
  #         legend.title = element_text(face = "bold")) +
  #   labs(x = "Max. future frequency >10%", 
  #        y = "Log10(no. libaries detected)")
  # 
  # merged %>%
  #   ggplot(aes(y = log10(global_n + 1), x = abs_charge)) +
  #   geom_bin2d() +
  #   scale_fill_viridis_c(trans = "log10") +
  #   geom_smooth(method = "lm", color = "black", fill = "grey") +
  #   theme_bw() +
  #   theme(legend.position = "top",
  #         text = element_text(family = "sans"),
  #         axis.title = element_text(face = "bold"),
  #         legend.title = element_text(face = "bold")) +
  #   labs(x = "Abs. charge change", y = "Log10(no. GISAID sequences)",
  #        fill = "No. mutations") 
  # 
  # merged %>%
  #   ggplot(aes(y = log10(global_n + 1), x = abs_mw)) +
  #   geom_bin2d() +
  #   scale_fill_viridis_c(trans = "log10") +
  #   geom_smooth(method = "lm", color = "black", fill = "grey") +
  #   theme_bw() +
  #   theme(legend.position = "top",
  #         text = element_text(family = "sans"),
  #         axis.title = element_text(face = "bold"),
  #         legend.title = element_text(face = "bold")) +
  #   labs(x = "Abs. charge change", y = "Log10(no. GISAID sequences)",
  #        fill = "No. mutations") 
  # 
  # ggarrange(p1, p2, p3, p4, nrow = 1, align = "hv")
  # 
  # ggsave(str_glue("results/modelling_out/predictors.{dataset}.pdf"), 
  #        width = 14, height = 4)
  
#   # Linear model
#   linreg <- lm(log10(global_n + 1) ~ n + max_freq + blosum62_score + abs_mw + abs_charge + abs_hydropathy, 
#                data = merged)
#   
#   res <- summary(linreg)
#   
#   prop_var <- relaimpo::calc.relimp(linreg)$lmg * 100
#   lin_res <- as.data.frame(coefficients(summary(linreg))) %>%
#     rownames_to_column("predictor") %>%
#     left_join(tibble(predictor = names(prop_var),
#                      perc_var_expl = prop_var) %>%
#                 arrange(desc(perc_var_expl))) %>%
#     mutate_if(is.numeric, function(x) {signif(x, 2)}) %>%
#     arrange(desc(perc_var_expl)) %>%
#     bind_rows(tibble(predictor = "total", perc_var_expl = res$adj.r.squared))
#   
#   fwrite(lin_res, str_glue("results/modelling_out/modelling_stats/linreg.{dataset}.csv"))
#   
#   # Logistic regression
#   logreg <- glm(future_high ~ log10(n) + max_freq + 
#                   blosum62_score + abs_mw + abs_charge + abs_hydropathy,
#                 data = merged,
#                 family = "binomial")
#   
#   logreg_res <- summary(logreg)
#   
#   parsed_logreg_res <- as.data.frame(coef(logreg_res)) %>%
#     mutate(OR = exp(Estimate)) %>%
#     mutate(across(everything(), function(x){signif(x, 3)}))
#   
#   fwrite(parsed_logreg_res, str_glue("results/modelling_out/modelling_stats/logreg.{dataset}.csv"))
#   
#   # Analysis of deviance
#   aod <- data.frame(anova(logreg))
#   null_deviance <- aod["NULL", "Resid..Dev"]
#   
#   aod %>%
#     dplyr::select(Deviance) %>%
#     mutate(dev_explained = Deviance / null_deviance * 100) 
#   
#   (null_deviance - logreg$deviance) / null_deviance * 100
#   
#   exp(coefficients(logreg))
# }

