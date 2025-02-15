rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")
gff3 <- fread("data/genomes/MN908947.3.gff3") %>%
  filter(V3 == "CDS") %>%
  separate(V9, c(NA, "gene_name"), "product\\=")

protein_hookup <- hookup %>%
  filter(codon_from_gene_start != -1) %>%
  distinct(protein_name, nucleotide_pos) %>%
  filter(protein_name != "NSP11")

file_dir <- "results/pipeline_out.020225/ivar_out/"
fna <- readDNAStringSet("data/alignments/reassembled.masked.filt.aln")

file_list <- list.files(file_dir, full.names = T)

# Check files
done <- gsub(".tsv", "", file_list)
done <- gsub(file_dir, "", done)
sum(names(fna) %in% done) == length(fna)

morsels <- foreach(file_name = file_list) %do% {
  id <- gsub(".tsv", "", file_name)
  id <- gsub(file_dir, "", id)
  
  fread(file_name) %>%
    mutate(id = id) %>%
    mutate(across(everything(), as.character))
}

# Parse frequency table
merged <- bind_rows(morsels) %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  mutate(pos = as.numeric(pos)) %>%
  select(id, pos, ref, alt, 
         ref_dp, ref_rv, ref_qual,
         alt_dp, alt_rv, alt_qual, 
         alt_freq, total_dp, pval, 
         pass, ref_codon, alt_codon,
         ref_aa, alt_aa, pos_aa, gff_feature) %>%
  filter(ref_aa != alt_aa) %>%
  filter(pass == "TRUE") %>%
  separate(gff_feature, c("protein_name"), "\\:", remove = F) %>%
  mutate(protein_name = gsub("ORF", "NS", protein_name))

orf1ab <- merged %>%
  filter(protein_name == "NS1ab")

non_orf1ab <- merged %>%
  filter(protein_name != "NS1ab")

parsed_orf1ab <- orf1ab %>%
  select(-protein_name) %>%
  left_join(hookup %>% 
              filter(protein_name != "NSP11") %>%
              distinct(pos = nucleotide_pos, protein_name)) 

parsed <- bind_rows(parsed_orf1ab, non_orf1ab) %>%
  mutate(mutation_name = str_glue("{protein_name}_{ref_aa}{pos_aa}{alt_aa}"))

parsed %>%
  group_by(id, mutation_name) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

  
fwrite(parsed, "results/allele_frequency_out/ivar_out/ivar_freq.all.raw.csv")

