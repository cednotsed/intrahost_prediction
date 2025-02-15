rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

fna <- readDNAStringSet("data/alignments/reassembled.gap_filtered.pango_filtered.aln")
file_dir <- "results/pipeline_out.020225/codon_out/"
file_list <- list.files(file_dir, full.names = T)

# Check files
done <- gsub(".csv", "", file_list)
done <- gsub(file_dir, "", done)
sum(names(fna) %in% done) == length(fna)

morsels <- foreach(file_name = file_list) %do% {
  id <- gsub(".csv", "", file_name)
  id <- gsub(file_dir, "", id)
  
  fread(file_name) %>%
    mutate(id = id) %>%
    mutate(across(everything(), as.character))
}

# Parse frequency table
merged <- bind_rows(morsels) %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  rename_all(~gsub("\\(|\\)", "", .x))

nt_meta <- merged %>%
  distinct(nt_position_gene) %>%
  separate(nt_position_gene, c("start"), "-", remove = F) 

# Parse codon frequency table
parsed_merged <- merged %>%
  left_join(nt_meta) %>%
  mutate(across(all_of(c("start", "nt_start_position", "nt_end_position", "coverage", "mutant_frequency")),
                as.numeric)) %>%
  mutate(codon_number = (nt_start_position - start) / 3 + 1) %>%
  select(id, ref_codon, var_codon = mutant_codon, protein_name = `#gene`, nt_position_gene,
         nt_start_position, nt_end_position, ref_AA = ref_aa, codon_number, var_AA = mutant_aa,
         coverage, freq = mutant_frequency, mutation_type = mutant_type) %>%
  # Merge NSP12
  mutate(codon_number = ifelse(protein_name == "NSP12_2", codon_number + 9, codon_number)) %>%
  # Parse protein names
  mutate(protein_name = ifelse(grepl("NSP12", protein_name), "NSP12", protein_name)) %>%
  mutate(protein_name = gsub("ORF", "NS", protein_name)) %>%
  mutate(protein_name = case_when(protein_name == "S" ~ "Spike",
                                  protein_name == "NSP3a" ~ "NSP3",
                                  protein_name %in% c("NSP3a", "NSP5a", "NSP15a", 
                                                      "NS3a", "NS8a") ~ gsub("a", "", protein_name),
                                  T ~ protein_name)) %>%
  mutate(mutation_name = ifelse(mutation_type == "NS", 
                                str_glue("{protein_name}_{ref_AA}{codon_number}{var_AA}"),
                                str_glue("{protein_name}_{ref_AA}{codon_number}{var_AA}"))) %>%
  mutate(freq = freq / 100) %>%
  filter(id %in% names(fna)) %>%
  mutate(n_mutations = str_count(var_codon, "[A-Z]")) %>%
  mutate(nt_start_position = nt_start_position + 1,
         nt_end_position = nt_end_position + 1) # Quasitools is zero-based, convert to one-based

fwrite(parsed_merged, "results/allele_frequency_out/codon_out/codon_freq.all.raw.csv.gz")

parsed_merged <- fread("results/allele_frequency_out/codon_out/codon_freq.all.raw.csv.gz")

# parsed_merged %>%
#   distinct(id, mutation_name, coverage) %>%
#   group_by(id, mutation_name) %>%
#   summarise(n = n_distinct(coverage)) %>%
#   arrange(desc(n))

# Sum frequencies for each mutation
parsed_agg <- parsed_merged %>%
  group_by(id, mutation_name, protein_name, ref_AA, codon_number, var_AA, mutation_type, nt_start_position, nt_end_position) %>%
  summarise(freq = sum(freq), 
            coverage = unique(coverage), 
            n_codons = n_distinct(var_codon)) %>%
  ungroup()

fwrite(parsed_agg, "results/allele_frequency_out/codon_out/codon_freq.all.summed.csv.gz")

# Filter variants
parsed_filt <- parsed_agg %>%
  filter(mutation_type == "NS") %>%
  filter(ref_AA != "*") %>%
  filter(var_AA != "*") %>%
  filter(coverage >= 100) %>%
  filter(freq > 0.03)

fwrite(parsed_filt, "results/allele_frequency_out/codon_out/missense_freq.filt.csv.gz")

# merged_agg <- merged_filt %>%
#   group_by(mutation_name, ref_AA, var_AA) %>%
#   summarise(n = n_distinct(id),
#             n_filt = n_distinct(id),
#             median_freq = median(freq),
#             max_freq = max(freq)) %>%
#   arrange(desc(n))
# fwrite(merged_agg, out_path2)
# 
# 
# 
#   
#          