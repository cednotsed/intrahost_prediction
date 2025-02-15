rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)

# Master df
dataset <- "ba1"

obs_df <- fread("results/linkage_out/observed_linkage/observed_Dprime.before_pirola.gt1000.csv")
intra_df <- fread(str_glue("results/linkage_out/intrahost_linkage.all/{dataset}.n5.with_zeroes.csv")) %>%
  mutate(mut1 = ifelse(mutation1 > mutation2, mutation1, mutation2),
         mut2 = ifelse(mutation1 > mutation2, mutation2, mutation1)) %>%
  dplyr::rename(intra_corr = corr)

hookup_df <- tibble(mutation_name = unique(c(obs_df$mut1, obs_df$mut2, intra_df$mut1, intra_df$mut2))) %>%
  separate(mutation_name, c("protein_name", "mut"), "_", remove = F) %>%
  mutate(codon_number = parse_number(mut)) %>%
  mutate(codon_name = str_glue("{protein_name}_{codon_number}")) %>%
  distinct(mutation_name, codon_name)

obs_parsed <- obs_df %>% 
  left_join(hookup_df %>% select(mut1 = mutation_name, codon1 = codon_name)) %>%
  left_join(hookup_df %>% select(mut2 = mutation_name, codon2 = codon_name)) %>%
  filter(codon1 != codon2)

intra_parsed <- intra_df %>% 
  left_join(hookup_df %>% select(mut1 = mutation_name, codon1 = codon_name)) %>%
  left_join(hookup_df %>% select(mut2 = mutation_name, codon2 = codon_name)) %>%
  filter(codon1 != codon2)

dnds_df <- fread("data/external_datasets/pond_2023.parsed.csv") %>%
  mutate(codon_name = str_glue("{protein_name}_{codon_number}")) %>%
  mutate(class = ifelse(class == "Invariable", "Purifying", class))

obs_plot_df <- obs_parsed %>%
  group_by(codon1, codon2) %>%
  summarise(max_Dprime = max(Dprime)) %>%
  left_join(dnds_df %>%
              select(codon1 = codon_name, class1 = class)) %>%
  left_join(dnds_df %>%
              select(codon2 = codon_name, class2 = class)) %>% 
  mutate(class_comb = case_when((class1 == "Purifying" & class2 == "Neutral")| 
                                  (class2 == "Purifying" & class1 == "Neutral") ~ "Purifying-Neutral",
                                (class1 == "Purifying" & class2 == "Diversifying")| 
                                  (class2 == "Purifying" & class1 == "Diversifying") ~ "Purifying-Diversifying",
                                (class1 == "Neutral" & class2 == "Diversifying")| 
                                  (class2 == "Neutral" & class1 == "Diversifying") ~ "Neutral-Diversifying",
                                (class1 == "Purifying" & class2 == "Purifying") ~ "Purifying-Purifying",
                                (class1 == "Neutral" & class2 == "Neutral") ~ "Neutral-Neutral",
                                (class1 == "Diversifying" & class2 == "Diversifying") ~ "Diversifying-Diversifying")) %>%
  filter(!is.na(class_comb)) %>%
  mutate(is_linked = max_Dprime > 0.9) %>%
  group_by(is_linked, class_comb) %>%
  summarise(n = n()) %>%
  group_by(is_linked) %>%
  mutate(prop = n / sum(n)) 

intra_plot_df <- intra_parsed %>%
  group_by(codon1, codon2) %>%
  summarise(max_corr = max(intra_corr)) %>%
  left_join(dnds_df %>%
              select(codon1 = codon_name, class1 = class)) %>%
  left_join(dnds_df %>%
              select(codon2 = codon_name, class2 = class)) %>% 
  mutate(class_comb = case_when((class1 == "Purifying" & class2 == "Neutral")| 
                                  (class2 == "Purifying" & class1 == "Neutral") ~ "Purifying-Neutral",
                                (class1 == "Purifying" & class2 == "Diversifying")| 
                                  (class2 == "Purifying" & class1 == "Diversifying") ~ "Purifying-Diversifying",
                                (class1 == "Neutral" & class2 == "Diversifying")| 
                                  (class2 == "Neutral" & class1 == "Diversifying") ~ "Neutral-Diversifying",
                                (class1 == "Purifying" & class2 == "Purifying") ~ "Purifying-Purifying",
                                (class1 == "Neutral" & class2 == "Neutral") ~ "Neutral-Neutral",
                                (class1 == "Diversifying" & class2 == "Diversifying") ~ "Diversifying-Diversifying")) %>%
  filter(!is.na(class_comb)) %>%
  mutate(is_linked = max_corr > 0.9) %>%
  group_by(is_linked, class_comb) %>%
  summarise(n = n()) %>%
  group_by(is_linked) %>%
  mutate(prop = n / sum(n),
         total = sum(n))

obs_plt <- obs_plot_df %>%
  ggplot(aes(x = class_comb, y = prop, fill = is_linked)) +
  geom_col(position = "dodge",
           color = "black") +
  scale_fill_manual(values = c("tan", "darkolivegreen")) +
  labs(x = "SAV pair classification", y = "Prop. pairs", fill = "Consensus linkage (D') >0.9") +
  theme(legend.position = "right") +
  ylim(0, 0.7) +
  theme_bw()
  

intra_plt <- intra_plot_df %>%
  ggplot(aes(x = class_comb, y = prop, fill = is_linked)) +
  geom_col(position = "dodge",
           color = "black") +
  scale_fill_manual(values = c("lightyellow1", "darkorchid4")) +
  labs(x = "SAV pair classification", y = "Prop. pairs", fill = "Intrahost linkage (r) > 0.9") +
  theme(legend.position = "right") +
  ylim(0, 0.7) +
  theme_bw()

ggarrange(obs_plt, intra_plt, nrow = 2)

intra_plot_df
obs_plot_df
