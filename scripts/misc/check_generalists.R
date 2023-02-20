# check abundance and prevalence of generalists in the myxo asv dataset

# load in packages
library(phyloseq)
library(tidyverse)

# read in phyloseq object
ps <- readRDS('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_asv_filt.rds')

# read in habitat preference
d_habpref <- read.csv('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_asv.csv')

# read in phyloseq object and grab tax table
d_taxa <- ps %>%
  phyloseq::tax_table() %>%
  data.frame() %>%
  janitor::clean_names() %>%
  rownames_to_column('otu')

# create d_meta
d_meta <- left_join(select(d_habpref, otu, habitat_preference, num_present), select(d_taxa, otu:family))

# filter for the rare traits
rare_traits <- filter(d_meta, str_detect(habitat_preference, 'marine_mud:terrestrial'))

# get data from phyloseq object
d_ps <- psmelt(ps) %>%
  janitor::clean_names()

d_ps <- group_by(d_ps, sample) %>%
  mutate(prop = abundance/sum(abundance)) %>%
  left_join(., select(d_meta, otu, habitat_preference))

filter(d_ps, otu %in% rare_traits$otu) %>% 
  filter(abundance > 0) %>%
  View()

d_summary <- filter(d_ps, abundance > 0) %>%
  group_by(otu, habitat_preference) %>%
  summarise(n = n(),
    across(prop, list(mean = mean, max = max, min = min)),
    .groups = 'drop')

# make plot looking at distribution of mean, max, and min proportions across all otus

d_summary %>%
  mutate(across(contains('prop'), log)) %>%
  pivot_longer(cols = n:prop_min, names_to = 'trait', values_to = 'value') %>%
  mutate(trait = gsub('abundance_', '', trait)) %>%
  ggplot(aes(habitat_preference, value)) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
  facet_wrap(~trait, scales = 'free') +
  theme_bw()

# basically yes although there are not many species that can live across all marine mud + terrestrial and freshwater + marine mud + terrestrial, they do not have lower average/min/max proportions than ASVs in other groups and are in at least 4 samples (which is the norm for lots of these myxococcus anyway).

# we could plot the outcome of their bootstraps to see what their distribution is like


  


