# Assign habitat preference to each OTU for each phyloseq object #

#--------------------------#
# what this script does ####
#--------------------------#

# reads in phyloseq object
# reads in distinct clusters assigned from the 16S amplicon data
# for each OTU it assigns a habitat preference based on its presence/absence in each cluster
# saves these things out

# load in packages 
library(phyloseq)
library(tidyverse)
library(speedyseq)
library(here)

# set seed to always get the same answer
set.seed(42)

# set where we are
here::i_am('scripts/sequencing_rpoB/processing/assign_habitat_preference_all.R')

# set percent similarity - those used in asvs_to_otus.R
percent_similarity <- c(99:90, 97.7, 85, 80, 'asv')

# read in group clusters from the 16S analysis
clusters <- read.csv('data/sequencing_16s/sample_cluster_assignments.csv') %>%
  mutate(across(where(is.numeric), as.character))

# we will use medoid clustering
clusters <- clusters %>%
  dplyr::select(., sample, medoid_nbclust) %>%
  mutate(clusters, clust = case_when(medoid_nbclust == '1' ~ 'terrestrial',
                                     medoid_nbclust == '2' ~ 'freshwater',
                                     medoid_nbclust == '3' ~ 'marine_mud'))

# from the PCoA plot of the myxobacteria, it can be seen that sample s46 actually is not terrestrial (as it was for the 16s)
# it is instead marine mud and the misidentification likely happened during the DNA extraction / sequencing

# change this cluster assignment
clusters <- mutate(clusters, clust = ifelse(sample == 'sample_s46', 'marine_mud', clust))

# look at availability of clusters
d_habitats <- clusters %>% 
  dplyr::select(sample, clust) %>%
  distinct() %>%
  group_by(clust) %>%
  tally() %>%
  mutate(prop_available = n / sum(n)) %>%
  dplyr::rename(num_available = n)

d_habitats
# terrestrial is more commonly sampled than mud and shore and freshwater is our least well sampled cluster

head(clusters)
clusters <- select(clusters, sample, clust)

#------------------------------#
# assign habitat preference ####
#------------------------------#

# rules and methods
# 1. If you are only present in a single cluster, you are assigned to that cluster
# 2. create dataframe for "availability" of clusters, essentially based on our sampling

# do this for each percent similarity object
for(i in 1:length(percent_similarity)){
  
  # define which level of OTU similarity we are using
  otu_similarity <- paste(percent_similarity[i], 'percent', sep = '')
  if(percent_similarity[i] == 'asv'){otu_similarity = 'asv'}
  
  # get number of raw taxa
  #raw_taxa_n <- readRDS(here(paste('data/sequencing_rpoB/phyloseq/myxococcus/clustered/ps_otu_', otu_similarity,  '.rds', sep = ''))) %>% ntaxa()
  
  # load in prevalence filtered myxococcus object
  ps_myxo <- readRDS(here(paste('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_', otu_similarity,  '_filt.rds', sep = '')))
  prev_taxa_n <- ntaxa(ps_myxo)
  
  # add in group clusters into this ps object #
  meta <- sample_data(ps_myxo) %>% data.frame() %>%
    mutate(sample = paste('sample_s', id, sep = '')) %>%
    left_join(., clusters)
  row.names(meta) <- paste('sample_s', meta$id, sep = '')
  sample_data(ps_myxo) <- sample_data(meta)
  
  # get dataset out of phyloseq
  d_ps <- psmelt(ps_myxo) %>%
    janitor::clean_names() %>%
    group_by(otu) %>%
    mutate(total_otu_abundance = sum(abundance),
           prop_of_all_otu = abundance / total_otu_abundance) %>%
    ungroup() %>%
    group_by(sample) %>%
    mutate(rel_abund = abundance / sum(abundance)) %>%
    ungroup()
  
  # calculate sample summary stats
  d_sample <- d_ps %>%
    filter(abundance > 0) %>%
    group_by(sample) %>%
    mutate(sample_depth = sum(abundance),
           n_species = n()) %>%
    ungroup() %>%
    select(sample, habitat_group_16s, location, clust, sample_depth, n_species) %>% 
    distinct()
  
  #------------------------------#
  # define habitat preference ####
  #------------------------------#
  
  # calculate statistics for each otu
  # 1. average abundance in each cluster when present
  # 2. number of times otu is present in each cluster
  # 3. average proportion in cluster when present
  # 4. proportion of presences in each cluster
  # calculate average abundance and average proportion
  d_pref <- filter(d_ps) %>%
    filter(abundance > 0) %>%
    group_by(otu, clust) %>%
    summarise(num_present = n(), .groups = 'drop',
              average_abundance = mean(abundance),
              average_prop = mean(rel_abund)) %>%
    group_by(otu) %>%
    mutate(habitats_present = n(),
           prop_present = num_present/sum(num_present)) %>%
    ungroup() %>%
    left_join(., d_habitats)
  
  
  # assign preference when an OTU is only present in one habitat cluster #
  
  # cannot do much more right now apart from assign to the cluster
  d_pref_single <- filter(d_pref, habitats_present == 1) %>%
    select(otu, clust, num_present, prop_available, average_prop, average_abundance) %>%
    distinct() %>%
    mutate(., habitat_preference = clust)

  # assign preference when an OTU is present in more than one habitat cluster #
  
  # filter for when present in > 1 sample
  d_pref_multiple <- filter(d_pref, habitats_present > 1)
  
  # 1. bootstrap a distribution of presence/absences for each otu
  # 2. each time pick a new sample with replacement from the real presences
  # 3. pick a new sample of 100 - this is arbitrary
  # 4. calculate proportion of presences in each habitat for each bootstrap
  
  # set number of boots
  n_boots=1000
  
  # resample from the original data (length of each new dataset is 100)
  # calculate proportion present in each habitat cluster in each new bootstrapped replicate dataset
  boots <- d_pref_multiple %>%
    select(., otu, clust, num_present) %>%
    uncount(weights = num_present) %>%
    group_by(otu) %>%
    dplyr::slice(rep(1:n(), times = n_boots)) %>%
    mutate(boot = rep(1:n_boots, each = n()/n_boots)) %>%
    group_by(otu, boot) %>%
    nest(data = clust) %>%
    mutate(new_data = purrr::map(data, ~sample(.x$clust, 100, replace = TRUE))) %>%
    unnest(new_data) %>%
    dplyr::rename(clust = new_data) %>%
    group_by(otu, boot, clust) %>%
    summarise(num_present = n(), .groups = 'drop') %>%
    complete(otu, boot, clust,
             fill = list(num_present = 0)) %>%
    group_by(boot, otu) %>%
    mutate(prop_present = num_present/sum(num_present)) %>%
    ungroup()
  
  boots <- left_join(boots, d_habitats)
  
  # calculate which quantile the available proportion is in of the distribution of used proportions
  d_quantiles <- group_by(boots, otu, clust) %>%
    summarise(quantile = sum(prop_present <= unique(prop_available))/n(), .groups = 'drop',
              quantile_threshold = quantile(prop_present, 0.975))
  
  # define affinity
  # you either can or cannot live there - affinity or no affinity
  # if prop available is <= 97.5% quantile of bootstrapped prop present, you have no affinity for that habitat
  # in other words, if after bootstrapping, only 2.5% of samples occur in that habitat at least the proportion expected by chance, they are deemed to have no affinity for that habitat
  #if prop available is <= 97.5% quantile of bootstrapped prop present, you have no affinity for that habitat
  # we define your ability to live somewhere as whats left over after we say you do not like to live somewhere
  d_quantiles <- mutate(d_quantiles, habitat_preference = ifelse(quantile <= 0.975, clust, 'no affinity'))
  
  d_quantiles2 <- filter(d_quantiles, habitat_preference != 'no affinity')
  
  d_quantile_summary <- select(d_quantiles2, otu, clust, quantile_threshold)
  
  #  calculate number of preferences for each otu
  d_pref_multiple2 <- group_by(d_quantiles2, otu) %>%
    tally() %>%
    dplyr::rename(number_habitats = n)
  
  # calculate summary of habitats
  boots_summary <- group_by(d_quantiles2, otu) %>%
    arrange(clust) %>%
    summarise(habitat_number = n(),
              habitat_preference = paste0(clust, collapse = ':'),
              .groups = 'drop')
  
  d_quantile_summary <- left_join(d_quantile_summary, boots_summary)
  
  # combine habitat preference datasets
  d_pref_multiple2 <- d_pref_multiple %>%
    left_join(., boots_summary) %>%
    select(., otu, num_present, average_prop, average_abundance, habitat_preference, clust)
  
  d_pref_all <- bind_rows(d_pref_single, d_pref_multiple2)
  
  d_pref_all_summary <- group_by(d_pref_all, otu, habitat_preference) %>%
    summarise(average_prop = mean(average_prop),
              average_abundance = mean(average_abundance),
              habitats_present = n(),
              num_present = sum(num_present),
              .groups = 'drop') %>%
    mutate(habitat_preference_n = str_count(habitat_preference, ':') +1)
  
  #------------------------#
  # save everything out ####
  #------------------------#
  
  # save out the bootstraps (useful for plotting)
  saveRDS(boots, paste('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/bootstraps/habpref_boots_', otu_similarity, '.rds', sep = ''))
  
  # save out habitat preferences
  write.csv(d_pref_all_summary, paste('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_', otu_similarity, '.csv', sep = ''), row.names = FALSE)
  
  # sample stats
  # d_sample is first
  d_ps2 <- filter(d_ps, abundance > 0) %>%
    select(sample, otu) %>%
    left_join(., select(d_pref_all_summary, otu, habitat_preference)) %>%
    group_by(sample, habitat_preference) %>%
    tally() %>%
    pivot_wider(names_from = 'habitat_preference', values_from = 'n') %>%
    mutate(across(everything(), function(x) replace_na(x, 0)))
  
  left_join(d_sample, d_ps2) %>%
    write.csv(., paste('data/sequencing_rpoB/phyloseq/myxococcus/sample_summary/sample_stats_', otu_similarity, '.csv', sep = ''), row.names = FALSE)

}






