rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(doParallel)

ref <- readDNAStringSet("data/genomes/MN908947.3.fna")
ref_vec <- rep("N", width(ref[1]))

file_dir <- "results/pipeline_out.tonkin/vcf_out.biallelic/"
file_list <- list.files(file_dir, full.names = T)

# Done
done <- list.files("results/pipeline_out.tonkin/consensus_out/")
done <- gsub(".fna", "", done)
done <- str_glue("results/pipeline_out.020225/vcf_out.biallelic/{done}.bcftools.biallelic.vcf.gz")
to_do <- file_list[!(file_list %in% done)]

# to_do <- file_list

cl <- makeCluster(16)
registerDoParallel(cl)

foreach(file_name = to_do, 
        .packages = c("tidyverse", "Biostrings", "foreach", "data.table")) %dopar% {
  # file_name <- to_do[1]
  # file_name = file_list[grepl("SAMD00268100", file_list)]
  print(id)
  id <- gsub(file_dir, "", file_name)
  id <- gsub(".bcftools.biallelic.vcf.gz", "", id)
  
  # Coverage file
  cov_df <- fread(str_glue("results/pipeline_out.tonkin/stats_out/{id}.coverage.txt")) %>%
    mutate(coverage = as.numeric(coverage))
  
  if(cov_df$coverage > 90 & nrow(cov_df) > 0) {
    temp <- fread(file_name)
    colnames(temp)[10] <- "format_string"
    
    temp_filt <- temp %>%
      mutate(POS = as.numeric(POS)) %>%
      mutate(QUAL = as.numeric(QUAL)) %>%
      separate(format_string, c(NA, NA, "total_depth", "dp"), "\\:") %>%
      separate(dp, c("ref_dp", "alt_dp"), "\\,") %>%
      mutate(total_depth = as.numeric(total_depth),
             ref_dp = as.numeric(ref_dp),
             alt_dp = as.numeric(alt_dp)) %>%
      filter(!is.na(QUAL)) %>%
      filter(total_depth >= 10) %>%
      filter(!grepl("INDEL", INFO)) %>%
      select(POS, REF, ALT, QUAL, total_depth, ref_dp, alt_dp) %>%
      mutate(ref_freq = ref_dp / total_depth, 
             alt_freq = alt_dp / total_depth) %>%
      group_by(POS) %>%
      filter(ref_freq > 0.5 | alt_freq > 0.5) %>%
      mutate(consensus = ifelse(ref_freq > alt_freq, REF, ALT)) %>%
      distinct(POS, consensus)
    
    counts <- temp_filt %>%
      group_by(POS) %>%
      summarise(n = n()) %>%
      filter(n > 1)
    
    if(nrow(counts) > 0) {
      print(str_glue("{id}"))
    } else {
    
      # Write consensus sequence
      consensus <- ref_vec
      consensus[temp_filt$POS] <- temp_filt$consensus
        
      # Print consensus stats
      perc <- sum(consensus == "N") / length(consensus) * 100
      print(str_glue("%Ns = {perc}"))
      
      # Write consensus
      dna <- DNAStringSet(paste0(consensus, collapse = ""))
      names(dna) <- id
      writeXStringSet(dna, str_glue("results/pipeline_out.tonkin/consensus_out/{id}.fna"))
    }
  }  
}


