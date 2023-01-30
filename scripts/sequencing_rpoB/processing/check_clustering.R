# Look at impact of clustering on the number of OTUs assigned at the family level
# Look at the impact of this pre and post prevalence filtering

#--------------------------#
# what this script does ####
#--------------------------#

# reads in myxo asv object
# reads in alignment and creates phyloseq object

# load in packages
library(phyloseq)
library(dada2)
library(DECIPHER)
library(here)
library(tidyverse)
library(speedyseq)

# set number of processors to use
num_processors <- 4

# set where we are
here::i_am('scripts/sequencing_rpoB/processing/asvs_to_otus.R')

# read in phyloseq object
ps <- readRDS(here('data/sequencing_rpoB/phyloseq/myxococcus/clustered/ps_otu_asv.rds'))

# get taxonomy table from ps object
tax_table <- tax_table(ps) %>%
  data.frame() %>%
  rownames_to_column(var = 'tip_label') %>%
  janitor::clean_names()

# prevalence filter the raw ASV data
ps_sub <- microViz::tax_filter(ps, min_prevalence = 4, min_total_abundance = 100)

# set percent similarity
percent_similarity <- c(99:90, 97.7, 85, 80)
cut_off = (100-percent_similarity) / 100

# run a for loop to calculate how many OTUs have been assigned to the family level

# create an empty data frame to populate the results
d_results <- tibble(percent_similarity = percent_similarity) %>%
  mutate(n_family_assigned = NA,
         n_total = NA)

for(i in 1:length(cut_off))
{
  # load in previously clustered datasets
  ps0_sub <- readRDS(here(paste('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_', percent_similarity[i], 'percent_filt.rds', sep = '')))
  
  # add these into tax table dataset
  d_taxa <- tax_table(ps) %>%
    data.frame() %>%
    rownames_to_column('otu') %>%
    janitor::clean_names() %>%
    left_join(ps0_cluster)
  

  d_results$n_family_assigned[i] <- tax_table(ps0_sub) %>% data.frame() %>% filter(!is.na(Family)) %>% nrow()
  d_results$n_total[i] <- ntaxa(ps0_sub)
}

d_results <- mutate(d_results, prop = n_family_assigned/n_total)

asv_filt_prop <- tax_table(ps_sub) %>% data.frame() %>% filter(!is.na(Family)) %>% nrow()/ntaxa(ps_sub)

#------------#
# not ran ####
#------------#

# load in cluster assignments
ps0_cluster <- readRDS(here(paste('data/sequencing_rpoB/phyloseq/myxococcus/clusters/ps_otu', percent_similarity[i], 'clusters.rds', sep = ''))) %>%
  rownames_to_column('otu')

# check proportion assigned to the same family within each cluster
d_tax_check <- group_by(d_taxa, cluster, family) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(cluster) %>%
  mutate(prop = n/sum(n)) %>%
  ungroup() %>%
  filter(!is.na(family)) %>%
  filter(prop >= 0.5)

d_results$n_family_assigned_filt_95[i] <- filter(ps0_cluster, otu %in% taxa_names(ps0_sub)) %>% filter(cluster %in% d_tax_check$cluster) %>% nrow()
d_results$n_family_assigned_95[i] <- filter(ps0_cluster, otu %in% taxa_names(ps0)) %>% filter(cluster %in% d_tax_check$cluster) %>% nrow()