rm(list = ls())
setwd("C:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)
require(rentrez)

df <- fread("data/metadata/sra_metadata/biosample_metadata_210324.full.new.tsv")

df1 <- df %>%
  filter(grepl("Severe acute respiratory syndrome coronavirus 2|CoV-2|cov2", taxonomy_name, ignore.case = T))

df2 <- df %>%
  filter(grepl("Severe acute respiratory syndrome coronavirus 2|CoV-2|cov2", biosample_title, ignore.case = T))

df3 <- df %>%
  filter(grepl("Severe acute respiratory syndrome coronavirus 2|CoV-2|cov2", isolate, ignore.case = T))

df4 <- df %>%
  filter(grepl("Severe acute respiratory syndrome coronavirus 2|CoV-2|cov2", strain, ignore.case = T))

df_filt <- bind_rows(df1, df2, df3, df4) %>%
  distinct()

dates_filt <- df_filt %>%
  distinct(collection_date) %>% 
  filter(!grepl("/", collection_date)) %>%
  separate(collection_date, c("Y", "M", "D"), "\\-", remove = F) %>%
  mutate(Y = as.numeric(Y),
         M = as.numeric(M),
         D = as.numeric(D)) %>%
  filter(!is.na(Y)) %>%
  filter(!is.na(M)) %>% 
  filter(!is.na(D)) %>% 
  filter(Y == 2021 & M == 2| 
           Y == 2021 & M == 6| 
           Y == 2021 & M == 12| 
           Y == 2022 & M == 6| 
           Y == 2023 & M == 2| 
           Y == 2023 & M == 12)

final <- df_filt %>%
  filter(collection_date %in% dates_filt$collection_date) %>%
  # Remove samples with "Plasmid" in metadata
  filter(!grepl("plasmid", strain)) %>%
  filter(host != "Plasmid") %>%
  # Remove lab hosts
  filter(lab_host == "") %>%
  # Remove samples with "16S" in metadata
  filter(!grepl("16s|wastewater|environmental", biosample_title)) %>% 
  filter(!grepl("wastewater", isolation_source)) %>% 
  filter(grepl("anterior|nasal|clinical|patient|fecal|human|naso", isolation_source)) %>%
  # Remove samples with missing accessions
  filter(biosample_accession != "") %>%
  # Get only human samples
  filter(host == "Homo sapiens") %>%
  left_join(dates_filt)

final %>% 
  group_by(Y, M) %>%
  summarise(n = n_distinct(biosample))

# Save biosample accessions to extract SRA run accessions using sra utils
final %>%
  dplyr::select(biosample_accession) %>%
  filter(!(biosample_accession %in% done$Sample)) %>%
  fwrite("data/metadata/sra_metadata/filtered_sra_accessions.transition.txt", eol = "\n", col.names = F)

# Read run accessions
run_df <- fread("data/metadata/sra_metadata/runinfo.transition.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  select(bioproject, biosample, biosample_accession = sample, 
         experiment, run, assemblyname,
         spots, bases, avglength, 
         platform, librarystrategy, libraryselection, 
         librarysource, librarylayout, download_path)

# Get missing efetch queries
final %>%
  left_join(run_df) %>%
  filter(is.na(librarysource)) %>%
  select(biosample_accession) %>%
  fwrite("data/metadata/sra_metadata/filtered_sra_accessions.transition.missing.txt", eol = "\n", col.names = F)
# fwrite("data/metadata/sra_metadata/filtered_sra_accessions.csv")

# Parse metadata
parsed <- final %>%
  left_join(run_df) %>%
  select(bioproject, biosample,
         experiment, run, assemblyname,
         collection_date, geo_loc_name, lat_lon,
         host, isolation_source, platform,
         spots, bases, avglength,
         librarystrategy, libraryselection, librarysource,
         librarylayout, download_path)

count_df <- parsed %>%
  group_by(biosample, bioproject) %>%
  summarise(n_platform = n_distinct(platform),
            n_layout = n_distinct(librarylayout),
            n_runs = n_distinct(run)) 

# Remove double layouts, double platform and double run samples
count_filt <- count_df %>%
  filter(n_platform == 1) %>%
  filter(n_layout == 1) %>%
  filter(n_runs == 1)

merged <- parsed %>%
  filter(biosample %in% count_filt$biosample) %>% 
  filter(platform == "ILLUMINA",
         librarylayout == "PAIRED") %>% 
  distinct(geo_loc_name, collection_date, .keep_all = T) %>%
  mutate(collection_month = format(as.Date(collection_date), "%Y-%m")) 

merged %>%
  fwrite("data/metadata/sra_metadata/filtered_sra_accessions.transition.csv")

merged %>%
  select(biosample) %>%
  fwrite("data/metadata/sra_metadata/filtered_sra_accessions.transition.accessions_only.paired.csv",
         eol = "\n",
         col.names = F)

