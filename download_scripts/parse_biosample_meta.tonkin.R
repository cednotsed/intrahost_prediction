rm(list = ls())
setwd("C:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)

df <- fread("data/metadata/sra_metadata/biosample_metadata_210324.full.new.tsv")
tonkin_bs <- fread("data/external_datasets/tonkin/tonkin.biosamples.txt", header = F)
runinfo <- fread("data/external_datasets/tonkin/runinfo.tonkin.csv") %>%
  rename_all(~tolower(gsub("\\ ", "_", .x)))

merged <- tibble(biosample = tonkin_bs$V1) %>%
  left_join(df)

test <- fread("data/metadata/all_sra_metadata.csv")
colnames(test)

parsed <- merged %>%
  rename_all(~tolower(gsub("\\ ", "_", .x))) %>%
  left_join(runinfo) %>%
  select(bioproject, biosample,
         experiment, run, assemblyname,
         collection_date, geo_loc_name, lat_lon,
         host, isolation_source, platform,
         spots, bases, avglength,
         librarystrategy, libraryselection, librarysource,
         librarylayout) %>% 
  separate(collection_date, c("Y", "M", "D"), "\\-", remove = F) %>%
  mutate(collection_month = ifelse(collection_date != "", str_glue("{Y}-{M}-01"), "2020-04-01")) %>%
  filter(librarylayout == "PAIRED") %>%
  mutate(sample_name = str_glue("{biosample}.{run}"))

fwrite(parsed, "data/external_datasets/tonkin/tonkin.metadata.csv")
