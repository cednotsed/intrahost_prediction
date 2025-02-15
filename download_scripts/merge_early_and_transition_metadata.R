rm(list = ls())
setwd("C:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)

early_meta <- fread("data/metadata/sra_metadata/filtered_sra_accessions.early.csv") %>%
  separate(collection_date, c("Y", "M"), "-", remove = F) %>%
  mutate(collection_month = str_glue("{Y}-{M}")) %>%
  dplyr::rename(platform = platforms) %>%
  filter(platform == "ILLUMINA",
         librarylayout == "PAIRED")

meta <- fread("data/metadata/sra_metadata/filtered_sra_accessions.transition.csv") %>%
  mutate(collection_date = as.character(collection_date)) %>%
  bind_rows(early_meta) %>%
  mutate(dataset = case_when(collection_month %in% c("2020-03", "2020-02", "2020-01", 
                                                     "2019-12", "2019-11", "2019-10") ~ "Early",
                             collection_month == "2021-02" ~ "Alpha",
                             collection_month == "2021-06" ~ "Delta",
                             collection_month == "2021-12" ~ "BA.1",
                             collection_month == "2022-06" ~ "BA.5",
                             collection_month == "2023-02" ~ "XBB",
                             collection_month == "2023-12" ~ "Pirola")) %>%
  mutate(collection_month = str_glue("{collection_month}-01")) %>%
  mutate(alias = tolower(gsub("\\.", "", dataset)))

meta %>%
  fwrite("data/metadata/all_sra_metadata.csv")

# meta %>%
#   select(biosample, collection_date) %>%
#   fwrite("data/metadata/all_sra_metadata.dates.tsv", sep = "\t")

