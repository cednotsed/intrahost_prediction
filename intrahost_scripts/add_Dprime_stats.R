rm(list = ls())
setwd("c:/git_repos/intrahost_prediction/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)

meta <- fread("data/metadata/all_sra_metadata.csv")
hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

datasets <- c("alpha", "delta", 
              "ba1", "ba5", "xbb",
              "pirola")

cl <- makeCluster(4)
registerDoParallel(cl)

morsels <- foreach(dataset = datasets, .packages = c("tidyverse", "foreach", "data.table")) %dopar% {
  print(dataset)
  # dataset = "delta"
  
  # Mutation stats
  stat_df <- fread(str_glue("results/mutation_stats/{dataset}.stats.csv"))
  
  # Prior fitness
  before_df <- fread(str_glue("results/allele_frequency_out/observed_out/before_{dataset}.aggregate.csv"))
  
  # Master df
  cor_df <- fread(str_glue("results/linkage_out/intrahost_linkage.all/{dataset}.n5.with_zeroes.csv")) %>%
    mutate(mut1 = ifelse(mutation1 > mutation2, mutation1, mutation2),
           mut2 = ifelse(mutation1 > mutation2, mutation2, mutation1)) %>%
    dplyr::rename(intra_corr = corr) %>%
    mutate(freq_range1 = mut1_max - mut1_min,
           freq_range2 = mut2_max - mut2_min)
  
  # Remove same codon mutations
  hookup_df <- tibble(mutation_name = unique(c(cor_df$mut1, cor_df$mut2))) %>%
    separate(mutation_name, c("protein_name", "mut"), "_", remove = F) %>%
    mutate(codon_number = parse_number(mut)) %>%
    mutate(codon_name = str_glue("{protein_name}_{codon_number}")) %>%
    select(mutation_name, codon_name)
  
  cor_parsed <- cor_df %>% 
    left_join(hookup_df %>% select(mut1 = mutation_name, codon1 = codon_name)) %>%
    left_join(hookup_df %>% select(mut2 = mutation_name, codon2 = codon_name)) %>%
    filter(codon1 != codon2) %>%
    arrange(desc(intra_corr))
  
  cor_filt <- cor_parsed %>%
    filter(mutation1 %in% stat_df$mutation_name) %>%
    filter(mutation2 %in% stat_df$mutation_name) %>%
    mutate(intra_corr = ifelse(intra_corr < 0, 0, intra_corr))
    # filter(freq_range1 > 0.5, freq_range2 > 0.5)
  
  cor_linked <- cor_filt %>%
    filter(intra_corr > 0.9)
  
  # Observed linkage stats
  obs_df <- fread(str_glue("results/linkage_out/observed_linkage/observed_linkage.before_{dataset}.gt1000.csv")) %>%
    mutate(mut1 = ifelse(mutation1 > mutation2, mutation1, mutation2),
           mut2 = ifelse(mutation1 > mutation2, mutation2, mutation1)) %>%
    dplyr::select(mut1, mut2, obs_corr = corr) %>%
    mutate(intra_corr = ifelse(obs_corr < 0, 0, obs_corr))
  
  obs_linked <- obs_df %>%
    filter(obs_corr > 0.9)
  
  dprime_df <- fread(str_glue("results/linkage_out/observed_linkage/observed_Dprime.before_{dataset}.gt1000.csv")) %>%
    mutate(Dprime = ifelse(Dprime < 0, 0, Dprime))
  
  dprime_linked <- dprime_df %>%
    filter(Dprime > 0.9)

  # Add max observed and intrahost correlation
  crumbs <- foreach(mut = unique(stat_df$mutation_name)) %do% {
    # mut = unique(stat_df$mutation_name)[1260]
    # mut <- "NS3_Q57H"
    # Intrahost
    temp <- cor_filt %>%
      filter(mutation1 == mut | mutation2 == mut)
    
    temp_linked <- cor_linked %>%
      filter(mutation1 == mut | mutation2 == mut)
    
    if(nrow(temp) == 0) {
      max_corr <- -1
      n_linked <- -1
      n_linked_binary <- -1
    } else {
      max_corr <- max(temp$intra_corr)
      n_linked <- nrow(temp_linked)
      n_linked_binary <- ifelse(nrow(temp_linked) > 0, 1, 0)
    }
    
    # Monthly freq correlation
    obs_temp <- obs_df %>%
      filter(mut1 == mut| mut2 == mut)
    
    obs_temp_linked <- obs_linked %>%
      filter(mut1 == mut| mut2 == mut)
    
    if(nrow(obs_temp) == 0) {
      max_obs_corr <- -1
      n_obs_linked <- -1
      n_obs_linked_binary <- -1
    } else {
      max_obs_corr <- max(obs_temp$obs_corr)
      n_obs_linked <- nrow(obs_temp_linked)
      n_obs_linked_binary <- ifelse(nrow(obs_temp_linked) > 0, 1, 0)
    }
    
    # Dprime
    temp_dprime <- dprime_df %>%
      filter(mut1 == mut | mut2 == mut)
    
    temp_dprime_linked <- dprime_linked %>%
      filter(mut1 == mut | mut2 == mut)
    
    linked_muts <- unique(c(temp_dprime_linked$mut1, temp_dprime_linked$mut2))
    linked_muts <- linked_muts[linked_muts != mut]
    fit_df <- tibble(mutation_name = linked_muts) %>%
      left_join(before_df) %>%
      mutate(global_n = replace_na(global_n, 0))
    
    if(nrow(temp_dprime) == 0) {
      max_dprime <- -1
      n_dprime_linked <- -1
      n_dprime_linked_binary <- -1
      median_dprime_fitness <- -1
      max_dprime_fitness <- -1
    } else {
      max_dprime <- max(temp_dprime$Dprime)
      n_dprime_linked <- nrow(temp_dprime_linked)
      n_dprime_linked_binary <- ifelse(nrow(temp_dprime_linked) > 0, 1, 0)
      median_dprime_fitness <- ifelse(nrow(temp_dprime_linked) > 0, median(fit_df$global_n), 0)
      max_dprime_fitness <- ifelse(nrow(temp_dprime_linked) > 0, max(fit_df$global_n), 0)
    }
    
    stat_df %>%
      filter(mutation_name == mut) %>%
      mutate(max_corr = max_corr, 
             n_linked = n_linked,
             n_linked_binary = n_linked_binary,
             max_obs_corr = max_obs_corr, 
             n_obs_linked = n_obs_linked,
             n_obs_linked_binary = n_obs_linked_binary,
             max_dprime = max_dprime,
             n_dprime_linked = n_dprime_linked,
             n_dprime_linked_binary = n_dprime_linked_binary,
             median_dprime_fitness = median_dprime_fitness,
             max_dprime_fitness = max_dprime_fitness)
  }
  
  link_stats <- bind_rows(crumbs)
  
  fwrite(link_stats, str_glue("results/mutation_stats/{dataset}.dprime_stats.csv"))
}

# link_stats %>%
#   filter(median_dprime_fitness != -1) %>%
#   arrange(log10(median_dprime_fitness)) %>%
#   ggplot(aes(x = log10(max_dprime_fitness), y = log10(global_n + 1))) +
#   geom_point() +
#   geom_smooth()
# 
# link_stats %>%
#   ggplot(aes(x = factor(n_linked), y = log10(global_n + 1))) +
#   geom_boxplot()
# 
# link_stats %>%
#   # filter(max_dprime != -1) %>%
#   ggplot(aes(x = max_dprime, y = log10(global_n + 1))) +
#   geom_bin2d() +
#   geom_smooth()
# # 
# dprime_df %>%
#   arrange(desc(Dprime)) %>% head() %>% View()
#   # filter(Dprime > 0) %>%
#   ggplot(aes(x = r2, y = Dprime ^ 2)) +
#   geom_bin2d()
#   
# dprime_df %>%
#   filter(Dprime != -1) %>%
#   summarise(neg = sum(Dprime < 0))
# 
# dprime_df <- fread(str_glue("results/linkage_out/observed_linkage/observed_Dprime.before_{dataset}.gt1000.csv"))
# dprime_df %>%
#   filter(Dprime != -1) %>%
#   summarise(neg = sum(Dprime < 0) / n())
# cor(dprime_df$r2, dprime_df$Dprime ^ 2, method = "spearman")
# link_stats %>% filter(global_n > 1e6) %>% View()