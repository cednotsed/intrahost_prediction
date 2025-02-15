rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)

data(BLOSUM62)

aa_df <- fread("data/metadata/aa_properties.csv")
possible_variants <- c("G", "A", "V", "L", "I",
                       "T", "S", "M", "C", "P",
                       "F", "Y", "W", "H", "K",
                       "R", "D", "E", "N", "Q")

changes <- as.data.frame(gtools::permutations(length(possible_variants), 2, possible_variants))
colnames(changes) <- c("ref_AA", "var_AA")

change_morsels <- foreach(i = seq(nrow(changes))) %do% {
  # i = 42
  print(i)
  row <- changes[i, ]
  
  row %>%
    left_join(aa_df %>% dplyr::select(ref_AA = amino_acid, ref_mw = mw, ref_charge = charge, ref_hydropathy = hydropathy)) %>%
    left_join(aa_df %>% dplyr::select(var_AA = amino_acid, var_mw = mw, var_charge = charge, var_hydropathy = hydropathy)) %>%
    mutate(delta_mw = var_mw - ref_mw,
           delta_charge = var_charge - ref_charge,
           delta_hydropathy = var_hydropathy - ref_hydropathy) %>%
    mutate(abs_mw = abs(delta_mw),
           abs_charge = abs(delta_charge),
           abs_hydropathy = abs(delta_hydropathy)) %>%
    mutate(blosum62_score = BLOSUM62[row$ref_AA, ifelse(row$var_AA == "stop", "*", row$var_AA)])
}

change_df <- bind_rows(change_morsels) %>%
  select(ref_AA, var_AA, blosum62_score,
         delta_charge, delta_mw, delta_hydropathy, 
         abs_charge, abs_mw, abs_hydropathy)

fwrite(change_df, "data/metadata/aa_properties.blosum62.csv")
