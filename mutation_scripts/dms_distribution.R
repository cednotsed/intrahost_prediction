rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(ggpubr)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

parsed_agg <- fread("results/mutation_out/monthly_freq_aggregate.csv")

possible_variants <- c("G", "A", "V", "L", "I",
                       "T", "S", "M", "C", "P",
                       "F", "Y", "W", "H", "K",
                       "R", "D", "E", "N", "Q")

# Add scores
aa_df <- fread("data/metadata/dms/final_variant_scores.csv") %>%
  filter(target == "Wuhan-Hu-1") %>%
  mutate(mutation_name = str_glue("Spike_{mutation}"))

abx_df <- fread("data/metadata/dms/MAP_paper_antibodies_raw_data.csv") %>%
  mutate(mutation_name = str_glue("Spike_{wildtype}{site}{mutation}")) %>%
  select(mutation_name, condition, mut_escape) %>%
  group_by(mutation_name) %>%
  summarise(mean_escape = sum(mut_escape) / n())
  # pivot_wider(id_cols = mutation_name, names_from = condition, values_from = mut_escape) %>%
  # column_to_rownames("mutation_name")

plot_df <- parsed_agg %>%
  mutate(is_fixed = max_prop > 0.9) %>%
  inner_join(aa_df %>% select(mutation_name, delta_bind, delta_expr)) %>%
  inner_join(abx_df)

p1 <- plot_df %>%
  ggplot(aes(x = delta_bind, y = log10(global_n))) +
  geom_bin2d() +
  scale_fill_viridis_c(trans = "log10") +
  geom_smooth(method = "lm", color = "black", fill = "grey") +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(family = "sans"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  labs(x = "Change in ACE2 binding", y = "Log10(no. GISAID sequences)",
       fill = "No. mutations") 

p2 <- plot_df %>%
  ggplot(aes(x = delta_expr, y = log10(global_n))) +
  geom_bin2d() +
  scale_fill_viridis_c(trans = "log10") +
  geom_smooth(method = "lm", color = "black", fill = "grey") +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(family = "sans"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  labs(x = "Change in RBD expression", y = "Log10(no. GISAID sequences)",
       fill = "No. mutations") 

ggarrange(p1, p2, nrow = 1)
ggsave("results/mutation_out/dms_v_fitness.png", width = 8, height = 4)

plot_df %>%
  ggplot(aes(x = delta_bind)) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "Hydropathy score", y = "Distinct mutations observed") +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        legend.position = "none",
        axis.title = element_text(face = "bold"))

plot_df %>%
  ggplot(aes(x = delta_expr)) +
  geom_histogram(fill = "slateblue", color = "black") +
  labs(x = "Charge", y = "Distinct mutations observed") +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        legend.position = "none",
        axis.title = element_text(face = "bold"))

plot_df %>%
  ggplot(aes(x = mean_escape)) +
  geom_histogram(fill = "darkseagreen", color = "black") +
  labs(x = "Molecular weight", y = "Distinct mutations observed") +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        legend.position = "none",
        axis.title = element_text(face = "bold"))

ggarrange(p1, p2, p3, nrow = 3, align = "hv")

ggsave("results/mutation_out/biochemistry.pdf", dpi = 600, height = 4, width = 4)

# Model
hydro <- MASS::rlm(delta_bind ~ is_fixed,
                   data = plot_df)

charge <- MASS::rlm(abs_charge ~ protein_name + is_fixed,
                    data = plot_df)

mw <- MASS::rlm(abs_mw ~ protein_name + is_fixed,
                data = plot_df)

sfsmisc::f.robftest(hydro, var = "is_fixedTRUE")
sfsmisc::f.robftest(charge, var = "is_fixedTRUE")
sfsmisc::f.robftest(mw, var = "is_fixedTRUE")

plot_df %>%
  ggplot(aes(x = is_fixed, y = abs_hydropathy)) +
  geom_boxplot() +
  geom_pwc()

plot_df %>%
  ggplot(aes(x = is_fixed, y = abs_charge)) +
  geom_boxplot() +
  geom_pwc()

plot_df %>%
  ggplot(aes(x = is_fixed, y = abs_mw)) +
  geom_boxplot() +
  geom_pwc()
# Model relationship to global_n
rlinreg <- MASS::rlm(log10(global_n) ~ abs_charge + abs_mw + abs_hydropathy,
                     data = plot_df,
                     psi = MASS::psi.bisquare)
summary(rlinreg)
sfsmisc::f.robftest(rlinreg, var = "abs_charge")
sfsmisc::f.robftest(rlinreg, var = "abs_mw")
sfsmisc::f.robftest(rlinreg, var = "abs_hydropathy")

cor.test(plot_df$global_n, plot_df$mean_escape, method = "spearman")
cor.test(plot_df$global_n, plot_df$delta_bind, method = "spearman")
cor.test(plot_df$global_n, plot_df$abs_hydropathy, method = "spearman")

plot_df %>%
  group_by(is_fixed) %>%
  summarise(median_charge = median(abs_charge),
            median_mw = median(abs_mw),
            median_hydropathy = median(abs_hydropathy))

# # Counts
# plot_df %>%
#   select(mutation_name, abs_mw, abs_charge, abs_hydropathy) %>%
#   pivot_longer(!mutation_name, names_to = "phenotype", values_to = "value") %>%
#   ggplot(aes(y = value, fill = phenotype)) +
#   geom_density() +
#   facet_wrap(~phenotype)
# 
# plot_df %>%
#   ggplot(aes(x = abs_mw)) +
#   geom_histogram(bins = 100)
# plot_df %>%
#   ggplot(aes(x = abs_charge)) +
#   geom_histogram(bins = 100)
# 
# plot_df %>%
#   filter(is_fixed) %>%
#   ggplot(aes(x = delta_hydropathy)) +
#   geom_density()




plot_df %>%
  ggplot(aes(x = abs_charge, fill = is_fixed)) +
  geom_histogram() +
  facet_grid(rows = vars(is_fixed),
             scale = "free")

plot_df %>%
  group_by(is_fixed) %>%
  summarise(median_charge = median(abs_charge),
            median_mw = median(abs_mw),
            median_hydropathy = median(abs_hydropathy))

plot_df %>%
  group_by(is_fixed) %>%
  summarise(median_charge = median(delta_charge),
            median_mw = median(delta_mw),
            median_hydropathy = median(delta_hydropathy))



summary(linreg)

plot_df %>%
  ggplot(aes(x = is_fixed, y = abs_hydropathy)) +
  geom_boxplot() +
  geom_pwc()
group_by(is_fixed) %>%
  summarise(median_charge = median(abs_mw))

plot_df %>%
  ggplot(aes(x = abs_hydropathy, y = log10(global_n))) +
  geom_bin2d()
facet_grid(rows = vars(is_fixed),
           scales = "free")
labs(x = "BLOSUM100 score", y = "Log10(no. genomes)") +
  scale_fill_manual(values = c("indianred", "turquoise3")) +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        legend.position = "none",
        axis.title = element_text(face = "bold"))

ggsave("results/mutation_out/global_n_versus_blosum100.pdf", dpi = 600, width = 4, height = 4)

linreg <- lm(log10(global_n) ~ protein_name + delta_charge + delta_mw + delta_hydropathy,
             data = plot_df)

rlinreg <- MASS::rlm(log10(global_n) ~ protein_name + delta_charge + delta_mw + delta_hydropathy,
                     data = plot_df,
                     psi = MASS::psi.bisquare)
sfsmisc::f.robftest(rlinreg, var = "delta_charge")

summary(rlinreg)
tibble(x = rlinreg$wresid) %>%
  ggplot(aes(x = x)) +
  geom_density()
coef(summary(linreg))

rlinreg <- glm(global_n ~ protein_name + abs_mw + abs_charge + abs_hydropathy,
               data = plot_df,
               family = "quasipoisson")

hist(linreg$residuals)
summary(linreg)

plot_df %>%
  mutate(is_fixed = max_prop > 0.9) %>%
  ggplot(aes(x = abs_charge, y = log10(global_n), color = is_fixed)) +
  geom_bin_2d() +
  geom_smooth(method = "lm")


plot_df %>%
  ggplot(aes(abs_charge)) +
  geom_histogram()
plot_df %>%
  group_by(blosum100_score) %>%
  summarise(n = n_distinct(mutation_name)) %>%
  ggplot(aes(x = factor(blosum100_score), y = n, fill = factor(blosum100_score))) +
  geom_bar(stat = "identity",
           color = "black") +
  scale_fill_viridis_d() +
  labs(x = "BLOSUM100 score", y = "Distinct mutations observed") +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        legend.position = "none",
        axis.title = element_text(face = "bold"))

ggsave("results/mutation_out/blosum_counts.pdf", dpi = 600, width = 4, height = 4)

# as.data.frame(blosum100) %>%
#   rownames_to_column("var1") %>%
#   pivot_longer(!var1, names_to = "var2", values_to = "score") %>%
#   filter(var1 != "")
#   arrange(desc(score)) %>% 

hookup %>%
  distinct(region, protein_name, codon_number) %>%
  filter(codon_number != -1) %>% 
  group_by(region) %>%
  summarise(n = n()) %>% View()
plot_df %>%
  mutate(ORF = ifelse(grepl("NSP", protein_name), "ORF1ab", protein_name)) %>%
  group_by(ORF) %>%
  summarise(n = n()) %>%
  View()

hookup %>% 
  filter(protein_name == "Spike") %>% View()

plot_df %>%
  mutate(conservative = blosum100_score >= 0) %>%
  group_by(is_fixed, conservative) %>%
  summarise(n = n()) %>%
  group_by(is_fixed) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = is_fixed, y = prop, fill = conservative)) +
  geom_bar(stat = "identity")


mat <- plot_df %>%
  mutate(conservative = blosum100_score >= 0) %>%
  group_by(is_fixed, conservative) %>%
  summarise(n = n()) %>%
  pivot_wider(id_cols = is_fixed, names_from = conservative, values_from = n) %>%
  column_to_rownames("is_fixed")

fisher.test(mat)

plot_df %>% 
  filter(!is.na(mean_escape)) %>%
  m
  filter(!is.na(delta_expr)) %>%
  summarise(prop = sum(delta_expr < 0) / n())
