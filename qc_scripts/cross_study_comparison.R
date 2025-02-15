rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")
fna <- readDNAStringSet("data/alignments/reassembled.masked.filt.aln")
meta <- fread("data/metadata/all_sra_metadata.csv")
sra_meta <- fread("data/metadata/all_sra_metadata.submission_dates.tsv")

parsed <- meta %>%
  filter(biosample %in% names(fna)) %>% 
  left_join(sra_meta %>% select(biosample = BioSample, submission_date = SubmissionDate)) %>% 
  mutate(submission_date = as.Date(submission_date),
         collection_date = as.Date(collection_date, "%Y-%m-%d")) %>%
  mutate(time_diff = difftime(submission_date, collection_date, units = "weeks")) %>%
  separate(geo_loc_name, c("country"), sep = "\\:\\ ", remove = F)

codon_morsels <- foreach(d = c("early", "alpha", "delta", "ba1", "ba5", "xbb", "pirola")) %do% {
  temp_df <- fread(str_glue("results/allele_frequency_out/codon_out/{d}.raw.filt.csv"))
}

snp_morsels <- foreach(d = c("early", "alpha", "delta", "ba1", "ba5", "xbb", "pirola")) %do% {
  temp_df <- fread(str_glue("results/allele_frequency_out/snp_out/{d}.merged_snp_frequencies.filt.csv"))
}


codon_df <- bind_rows(codon_morsels) %>%
  filter(freq < 0.5) %>%
  filter(id %in% names(fna))

snp_df <- bind_rows(snp_morsels) %>%
  filter(freq < 0.5) %>%
  filter(id %in% names(fna))

codon_counts <- codon_df %>%
  group_by(id) %>%
  summarise(n_codons = n()) %>%
  ungroup()

snp_counts <- snp_df %>%
  group_by(id) %>%
  summarise(n_snps = n()) %>% 
  ungroup()

combined <- tibble(id = names(fna)) %>%
  left_join(snp_counts) %>%
  left_join(codon_counts) %>%
  mutate(n_snps = replace_na(n_snps, 0),
         n_codons = replace_na(n_codons, 0))

combined %>%
  summarise(median_snps = median(n_snps),
            median_codons = median(n_codons),
            mean_snps = mean(n_snps),
            mean_codons = mean(n_codons))
cor.test(combined$n_snps, combined$n_codons, method = "spearman")  

combined %>%
ggplot(aes(x = n_snps, y = n_codons)) +
  geom_point()

codon_counts %>%
  ungroup() %>%
  summarise(mean = median(n_codons))

protein_meta <- hookup %>%
  distinct(protein_name, protein_length) %>%
  filter(!is.na(protein_length)) %>%
  group_by(protein_name) %>%
  summarise(protein_length = sum(protein_length))

order_df <- hookup %>%
  filter(ref_AA != -1) %>%
  distinct(protein_name)

parsed <- codon_df %>%
  group_by(protein_name) %>%
  summarise(total_variants = n_distinct(mutation_name)) %>%
  left_join(protein_meta) %>%
  filter(protein_name != "NSP11") %>%
  mutate(ratio = total_variants / protein_length) %>%
  mutate(ratio = ratio / length(fna)) %>%
  mutate(protein_name = factor(protein_name, order_df$protein_name)) %>%
  arrange(desc(ratio))

parsed %>%
  arrange(desc(ratio)) %>%
  ggplot(aes(x = protein_name, y = ratio)) +
  geom_col()
   
# protein_meta
# snp_counts <- snp_df %>% 
#   group_by(id) %>%
#   summarise(n_snps = n())
# 
# snp_counts %>%
#   left_join(codon_counts) %>%
#   ggplot(aes(x = n_snps, y = n_codons)) +
#   geom_point()


# Lythgoe et al.
lythgoe_df <- fread("data/external_datasets/lythgoe_mutation_table.csv")

parsed_lythgoe <- lythgoe_df %>%
  mutate(protein_name = gsub("ORF", "NS", protein_name)) %>%
  mutate(protein_name = gsub("nsp", "NSP", protein_name)) %>%
  mutate(protein_name = gsub("\\*", "", protein_name)) %>%
  mutate(protein_name = gsub("NSP5A", "NSP5", protein_name)) %>%
  mutate(protein_name = gsub("NS3a", "NS3", protein_name)) %>%
  mutate(protein_name = ifelse(protein_name == "S", "Spike", protein_name)) %>%
  filter(protein_name != "NS1b") %>%
  mutate(protein_name = factor(protein_name, order_df$protein_name)) %>%
  mutate(ratio = non_syn / length) %>%
  mutate(ratio = ratio / 1390) %>%
  arrange(desc(ratio))

# Tonkin-Hill et al.
tonkin <- fread("data/external_datasets/tonkin-hill_table_s2.csv")
n_tonkin <- n_distinct(tonkin$sampleID)

protein_hookup <- hookup %>%
  filter(codon_from_gene_start != -1) %>%
  distinct(protein_name, nucleotide_pos)

parsed_tonkin <- tonkin %>%
  filter(vaf < 0.5) %>%
  filter(vaf > 0.03) %>%
  filter(impact == "Missense") %>% 
  left_join(protein_hookup %>% select(protein_name, pos = nucleotide_pos)) %>%
  group_by(protein_name) %>% 
  summarise(n_muts = n_distinct(aachange)) %>% 
  left_join(protein_meta) %>%
  mutate(ratio = n_muts / protein_length) %>%
  mutate(ratio = ratio / 1181)
  

bind_rows(parsed_lythgoe %>%
            mutate(study = "Lythgoe et al. 2021"), 
          parsed %>% 
            mutate(study = "Our study"),
          parsed_tonkin %>%
            mutate(study = "Tonkin-Hill et al. 2021")) %>%
  mutate(study = factor(study, c("Our study", "Tonkin-Hill et al. 2021", "Lythgoe et al. 2021"))) %>%
  mutate(protein_name = factor(protein_name, order_df$protein_name)) %>%
  ggplot(aes(x = protein_name, y =ratio)) +
  geom_col(color = "black") +
  labs(x = "Protein name", y = "SAVs / gene length") +
  facet_grid(rows = vars(study), scales = "free")

ggsave("results/qc_out/cross_study_comparison.pdf", dpi = 600, width = 12, height = 5)

parsed_tonkin %>%
  group_by(sampleID) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) 
