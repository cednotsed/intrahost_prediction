rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

bed <- fread("data/metadata/wuhan-hu-1_genome_annotations_V2.bed")


gff <- bed %>%
  mutate(type = "RefSeq", region = "CDS", score = ".", strand = "+", phase = ".") %>%
  mutate(attributes = str_glue("gene={V4}")) %>%
  select(V1, type, region, V2, V3, score, strand, phase, attributes)

fwrite(gff, "data/genomes/MN908947.3.parsed.gff3",
       eol = "\n",
       col.names = F,
       sep = "\t")  
