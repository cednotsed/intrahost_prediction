rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(randomcoloR)


df <- fread("data/external_datasets/pond_2023.csv") %>%
  rename_all(~tolower(gsub("\\ ", "_", .x)))

parsed <- df %>%
  separate(codon, c("protein_name", "codon_number"), "\\/") %>%
  filter(protein_name != "ORF5") %>%
  mutate(protein_name = case_when(protein_name == "leader" ~ "nsp1",
                                  protein_name == "3C" ~ "nsp5",
                                  protein_name == "RdRp" ~ "nsp12", 
                                  protein_name == "helicase" ~ "nsp13", 
                                  protein_name == "exonuclease" ~ "nsp14", 
                                  protein_name == "endornase" ~ "nsp15", 
                                  protein_name == "methyltransferase" ~ "nsp16",
                                  protein_name == "ORF3a" ~ "NS3",
                                  protein_name == "S" ~ "Spike",
                                  protein_name == "ORF3a" ~ "NS3",
                                  protein_name == "ORF6" ~ "NS6",
                                  protein_name == "ORF7a" ~ "NS7a",
                                  protein_name == "ORF8" ~ "NS8",
                                  protein_name == "ORF10" ~ "NS10",
                                  TRUE ~ protein_name)) %>%
  mutate(protein_name = gsub("nsp", "NSP", protein_name)) %>%
  dplyr::rename(nt_start = location)

parsed %>%
  fwrite("data/external_datasets/pond_2023.parsed.csv")

