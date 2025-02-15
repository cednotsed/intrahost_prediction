rm(list = ls())
setwd("C:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(combinat)

annot_df <- fread("data/metadata/wuhan-hu-1_genome_annotations_V2.csv") %>%
  filter(region_type != "non-coding")
  # mutate(protein_name = ifelse(grepl("NSP12", protein_name), "NSP12", protein_name))

fna <- readDNAStringSet("data/genomes/MN908947.3.fna")

proteins <- annot_df %>% distinct(protein_name)
proteins

# Get codon changes
codons <- names(GENETIC_CODE)
pairs <- permutations(length(codons), 2, codons, repeats.allowed = F)

# Filter only non-synonymous changes
change_morsels <- foreach(i = seq(nrow(pairs))) %do% {
  print(i)
  # i = 1000
  pair <- pairs[i, ]
  c1 <- pair[1]
  c2 <- pair[2]
  
  if (GENETIC_CODE[[c1]] != GENETIC_CODE[[c2]]) {
    # Compare characters
    seq1 <- str_split(c1, "")[[1]]
    seq2 <- str_split(c2, "")[[1]]
    change_temp <- as.data.frame(cbind(seq1, seq2)) %>%
      mutate(codon_position = seq(3)) %>%
      filter(seq1 != seq2) %>%
      mutate(codon = c1) %>%
      mutate(change = str_glue("{seq1}{codon_position}{seq2}"))
    
    tibble(ref_codon = c1, 
           var_codon = c2,
           changes = paste0(change_temp$change, collapse = "_"),
           n_changes = nrow(change_temp),
           ref_AA = GENETIC_CODE[[c1]], 
           var_AA = GENETIC_CODE[[c2]])
  }
}

change_df <- bind_rows(change_morsels) %>%
  # filter(n_changes > 1) %>%
  separate(changes, c("change1", "change2", "change3"), "_", remove = F) %>%
  mutate(cpos1 = parse_number(change1),
         cpos2 = parse_number(change2),
         cpos3 = parse_number(change3)) %>%
  mutate(change1 = gsub('[0-9]+', '', change1),
         change2 = gsub('[0-9]+', '', change2),
         change3 = gsub('[0-9]+', '', change3))

all_cds <- foreach(protein = proteins$protein_name) %do% {
  protein_annot <- annot_df %>%
    filter(protein_name == protein)
  
  cds <- substr(fna, protein_annot$genome_start, protein_annot$genome_end)
  protein_length <- nchar(cds) / 3
  codon_list <- sapply(seq(from = 1, to = nchar(cds), by = 3), function(i) substr(cds, i, i + 2))
  
  codon_df <- tibble(ref_codon = codon_list,
                     codon_number = seq(protein_length),
                     first_nuc_pos = seq(protein_annot$genome_start, 
                                         protein_annot$genome_end)[rep(c(T, F, F), 
                                                                       length(codon_list))]) %>%
    left_join(change_df) %>% 
    mutate(across(contains("cpos"), function(x){x + first_nuc_pos - 1})) %>%
    mutate(change1 = paste0(substr(change1, 1, 1), cpos1, substr(change1, 2, 2)),
           change2 = paste0(substr(change2, 1, 1), cpos2, substr(change2, 2, 2)),
           change3 = paste0(substr(change3, 1, 1), cpos3, substr(change3, 2, 2))) %>% 
    mutate(allele_changes = case_when(n_changes == 1 ~ change1,
                                      n_changes == 2 ~ str_glue("{change1}_{change2}"),
                                      n_changes == 3 ~ str_glue("{change1}_{change2}_{change3}"))) %>%
    mutate(protein_name = protein)
  
  return(codon_df)
}

final <- bind_rows(all_cds) %>%
  mutate(mutation_name = str_glue("{protein_name}_{ref_AA}{codon_number}{var_AA}"))

final %>%
  select(mutation_name, ref_AA, codon_number, 
         var_AA, ref_codon, var_codon, 
         changes, allele_changes, n_changes)
