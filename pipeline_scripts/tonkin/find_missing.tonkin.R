rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)

file_list <- list.files("results/pipeline_out.tonkin/stats_out", full.names = T)
file_list <- file_list[grepl("coverage", file_list)]
head(file_list)

done <- foreach(file_name = file_list, .combine = "c") %do% {
    print(file_name)
    temp <- fread(file_name)
    id <- gsub("results/pipeline_out.tonkin/stats_out/|.coverage.txt", "", file_name)

    if(nrow(temp) > 0) {
        #print(id)
        return(id)
    } else {
        #print("NO")
        return(NULL)
    }
}

merged <- tibble(id = done) %>%
  separate(id, c("biosample", "run")) %>%
  group_by(biosample) %>%
  summarise(n = n_distinct(run))

# parsed %>%
#      fwrite("data/metadata/sra_metadata/filtered_sra_accessions.transition.accessions_only.paired.missing.csv",
#             eol = "\n",
#             col.names = F)


