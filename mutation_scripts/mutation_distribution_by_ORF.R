rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(randomcoloR)
require(paletteer)

parsed_agg <- fread("results/allele_frequency_out/observed_out/all_time.aggregate.csv")

# Get protein lengths
protein_meta <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv") %>%
  filter(!is.na(protein_length)) %>%
  distinct(region, protein_length) %>%
  mutate(region = case_when(region == "ORF1ab" ~ "ORF1ab",
                            region == "S" ~ "Spike",
                            grepl("ORF", region) & region != "ORF1ab" ~ gsub("ORF", "NS", region),
                            TRUE ~ region)) %>%
  mutate(region = ifelse(region == "NS3a", "NS3", region))

orf1ab <- protein_meta %>%
  filter(region == "ORF1ab") %>%
  group_by(region) %>%
  summarise(region_length = sum(protein_length)) %>%
  mutate(region_length = 7097)

protein_parsed <- protein_meta %>%
  filter(region != "ORF1ab") %>%
  dplyr::rename(region_length = protein_length) %>%
  bind_rows(orf1ab) %>%
  # Remove start and stop codons
  mutate(region_length = region_length - 2) %>%
  mutate(possible_mutations = 20 * region_length)

fixed <- parsed_agg %>%
  filter(max_prop > 0.9)

drift <- parsed_agg %>%
  filter(max_prop < 0.1)

major <- parsed_agg %>%
  filter(max_prop > 0.1 & max_prop < 0.9)

nrow(fixed) / nrow(parsed_agg) * 100
nrow(drift) / nrow(parsed_agg) * 100
nrow(major) / nrow(parsed_agg) * 100

plot_df <- parsed_agg %>%
  left_join(protein_parsed) %>%
  mutate(mut_cat = case_when(mutation_name %in% fixed$mutation_name ~ ">90%",
                             mutation_name %in% drift$mutation_name ~ "<10%",
                             mutation_name %in% major$mutation_name ~ "10-90%")) %>%
  mutate(region = factor(region, c("ORF1ab", "Spike", "NS3",
                                   "E", "M", "NS6", 
                                   "NS7a", "NS7b", "NS8",
                                   "N", "NS10"))) %>%
  group_by(region, mut_cat, region_length) %>%
  summarise(n = n()) %>%
  mutate(norm_n = n / region_length)

pal <- c("#DCAC51", "#CB9E88", "#B44CE7", "#DEE561", "#DDAAD2", "#D5DCD9", "#D1E0A2", "#897FD4", "#84B2D8",
         "#7ADD91", "#DC6AC2", "#90EA4F", "#E45B63", "#70D7D0")

pal <- c("#5D74A5FF", "#B0CBE7FF", "#FEF7C7FF", 
         "#EBA07EFF", "#A8554EFF","#6C568CFF", 
         "#9386A6FF", "#BFCDD9FF", "#7F8C72FF",
         "#607345FF")

# pal <- setNames(pal, unique(plot_df$region_name))

plot_df %>%
  mutate(mut_cat = factor(mut_cat, c("<10%", "10-90%", ">90%"))) %>%
  ggplot(aes(x = region, y = norm_n, fill = region)) +
  geom_bar(stat = "identity",
           color = "black") +
  geom_text(aes(label = n),
            color = "grey90",
            angle = 90,
            position = position_stack(vjust = 0.5)) +
  facet_grid(rows = vars(mut_cat),
             scales = "free") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        legend.position = "none",
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  labs(x = "Region", y = "Distinct mutations / protein length")
 
ggsave("results/mutation_out/mutation_by_ORF.pdf", dpi = 600, width = 6, height = 4)

