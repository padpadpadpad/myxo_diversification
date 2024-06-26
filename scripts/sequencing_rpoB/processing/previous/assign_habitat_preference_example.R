# Assign habitat preference to each OTU for each phyloseq object #

#--------------------------#
# what this script does ####
#--------------------------#

# reads in phyloseq object
# reads in distinct clusters assigned from the 16S amplicon data
# for each OTU it assigns a habitat preference based on its presence/absence in each cluster
# saves these things out

# load in necessary packages
library(phyloseq)
library(tidyverse)
library(speedyseq)

#-------------------------------#
# extra functions to load in ####
#-------------------------------#

# function to sample n groups
sample_n_groups = function(tbl, size, replace = FALSE, weight = NULL) {
# regroup when done
grps = tbl %>% groups %>% lapply(as.character) %>% unlist
# check length of groups non-zero
keep = tbl %>% summarise() %>% ungroup() %>% sample_n(size, replace, weight)
# keep only selected groups, regroup because joins change count.
# regrouping may be unnecessary but joins do something funky to grouping variable
tbl %>% right_join(keep, by=grps) %>% group_by(.dots = grps)
}

# read in group clusters from the 16S analysis




#----------------------------#
# read in phyloseq object ####
#----------------------------#

# define which level of OTU similarity we are using
otu_similarity <- '97.7percent'

# get number of raw taxa
raw_taxa_n <- readRDS(paste('sequencing_rpoB/data/output/run_myxo_gtdbr202/ps_otu_', otu_similarity,  '.rds', sep = '')) %>% ntaxa()
  
# load in prevalence filtered myxococcus object
ps_myxo <- readRDS(paste('sequencing_rpoB/data/output/run_myxo_gtdbr202/prevalence_filtered/ps_myxo_', otu_similarity,  '.rds', sep = ''))
prev_taxa_n <- ntaxa(ps_myxo)

#----------------------------------------------#
# add in group clusters into this ps object ####
#----------------------------------------------#

clusters <- read.csv('sequencing_16S/data/analysis/asv_cluster_assignments.csv') %>%
  mutate(across(where(is.numeric), as.character))

head(clusters)

meta <- sample_data(ps_myxo) %>% data.frame() %>%
  mutate(sample = paste('sample_s', id, sep = '')) %>%
  left_join(., clusters)
row.names(meta) <- paste('sample_s', meta$id, sep = '')
sample_data(ps_myxo) <- sample_data(meta)

sample_data(ps_myxo)

# check number of reads in each sample ####
sample_sums(ps_myxo) %>% sort()
ps_myxo

#-----------------------------------------------#
# calculate summary data for phyloseq object ####
#-----------------------------------------------#

# three of the clustering methods. Hierarchical NB, hierarchical gap, and medoid clustering NB give almost exactly the same results. 3 clusters, broadly defined as terrestrial, freshwater, and mud and shore.

# get dataset out of phyloseq
d_ps <- psmelt(ps_myxo) %>%
  janitor::clean_names() %>%
  group_by(otu) %>%
  mutate(total_otu_abundance = sum(abundance),
         prop_of_all_otu = abundance / total_otu_abundance) %>%
  ungroup() %>%
  mutate(clust = case_when(medoid_nbclust == '2' ~ 'terrestrial',
                           medoid_nbclust == '3' ~ 'freshwater',
                           medoid_nbclust == '1' ~ 'mud_and_shore')) %>%
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
  select(sample, habitat_group, location, clust, sample_depth, n_species) %>% 
  distinct()

#------------------------------#
# define habitat preference ####
#------------------------------#

# rules and methods
# 1. If you are only present in a single cluster, you are assigned to that cluster
# 2. create dataframe for "availability" of clusters, essentially based on our sampling

# look at availability of clusters
d_habitats <- d_ps %>% 
  select(sample, clust) %>%
  distinct() %>%
  group_by(clust) %>%
  tally() %>%
  mutate(prop_available = n / sum(n)) %>%
  rename(num_available = n)

d_habitats
# terrestrial is more commonly sampled than mud and shore and freshwater is our least well sampled cluster

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

#-------------------------------------------------------------------------#
# assign preference when an OTU is only present in one habitat cluster ####
#-------------------------------------------------------------------------#

# cannot do much more right now apart from assign to the cluster
d_pref_single <- filter(d_pref, habitats_present == 1) %>%
  select(otu, clust, num_present, prop_available, average_prop, average_abundance) %>%
  distinct() %>%
  mutate(., habitat_preference = clust)
  
group_by(d_pref_single, habitat_preference) %>%
  tally() %>%
  mutate(prop = n/sum(n))
# number of single-use ASVs pretty close to that expected by availability

ggplot(d_pref_single, aes(num_present)) +
  facet_wrap(~clust) +
  geom_bar(col = 'black', fill = 'white') +
  theme_bw()

ggplot(d_pref_single, aes(average_prop)) +
  facet_wrap(~clust) +
  geom_histogram(col = 'black', fill = 'white') +
  theme_bw()
# most single habitat OTUs are only present in the bare minimum number of samples (4)

#------------------------------------------------------------------------------#
# assign preference when an OTU is present in more than one habitat cluster ####
#------------------------------------------------------------------------------#

# filter for when present in > 1 sample
d_pref_multiple <- filter(d_pref, habitats_present > 1)

# 1. bootstrap a distribution of presence/absences for each otu
# 2. each time pick a new sample with replacement from the real presences
# 3. pick a new sample of 100 - this is arbitrary
# 4. calculate proportion of presences in each habitat for each bootstrap

# set number of boots
n_boots=1000

# resample from the original data (length of each new dataset is 100)
boots <- d_pref_multiple %>%
  select(., otu, clust, num_present) %>%
  uncount(weights = num_present) %>%
  group_by(otu) %>%
  slice(rep(1:n(), times = n_boots)) %>%
  mutate(boot = rep(1:n_boots, each = n()/n_boots)) %>%
  group_by(otu, boot) %>%
  nest(data = clust) %>%
  mutate(new_data = map(data, ~sample(.x$clust, 100, replace = TRUE))) %>%
  unnest(new_data) %>%
  rename(clust = new_data) %>%
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
            quantile_threshold = quantile(prop_present, 0.975),
            prop_present_mean = mean(prop_present),
            prop_present_max = max(prop_present),
            prop_present_min = min(prop_present),
            prop_available = unique(prop_available))

# define affinity
# you either can or cannot live there - affinity or no affinity
# if prop available is <= 97.5% quantile of bootstrapped prop present, you have no affinity for that habitat
# in other words, if after bootstrapping, only 2.5% of samples occur in that habitat at least the proportion expected by chance, they are deemed to have no affinity for that habitat
#if prop available is <= 97.5% quantile of bootstrapped prop present, you have no affinity for that habitat
# we define your ability to live somewhere as whats left over after we say you do not like to live somewhere
d_quantiles <- mutate(d_quantiles, habitat_preference = ifelse(quantile <= 0.975, clust, 'no affinity'))

d_quantiles2 <- filter(d_quantiles, habitat_preference != 'no affinity')

unique(d_quantiles$otu) %>% length()
unique(d_quantiles2$otu) %>% length() # not lost any this is good

d_quantile_summary <- select(d_quantiles2, otu, clust, quantile_threshold)

#  calculate number of preferences for each otu
d_pref_multiple2 <- group_by(d_quantiles2, otu) %>%
  tally() %>%
  rename(number_habitats = n)

quantiles_to_plot <- select(d_quantiles2, otu, clust, quantile_threshold) %>%
  left_join(., d_pref_multiple2)
filter(quantiles_to_plot, number_habitats == 3)

boots2 <- left_join(boots, quantiles_to_plot)

# lets visualise all of our "generalists"
filter(boots2, number_habitats ==3) %>%
  ggplot(., aes(prop_present)) +
  facet_grid(otu ~clust) +
  geom_histogram(col = 'black', fill = 'white', bins = 30) +
  geom_point(aes(quantile_threshold,10), shape = 21, fill = 'grey', size = 3) +
  theme_bw() +
  geom_vline(aes(xintercept = prop_available), d_habitats, col = 'red')

# calculate summary of habitats
boots_summary <- group_by(d_quantiles2, otu) %>%
  arrange(clust) %>%
  summarise(habitat_number = n(),
    habitat_preference = paste0(clust, collapse = ':'),
    .groups = 'drop')

boots_summary %>% group_by(habitat_preference, habitat_number) %>%
  tally() %>%
  arrange(habitat_number)
# lots of freshwater:mud and shore and freshwater:terrestrial
# very few true generalists or mud and shore:terrestrial

quantiles_to_plot <- select(d_quantiles, otu, clust, quantile_threshold) %>%
  left_join(., boots_summary)

filter(quantiles_to_plot, otu == 'otu_5960')

boots2 <- left_join(boots, quantiles_to_plot)

filter(boots2, otu == 'otu_5960')

d_quantile_summary <- left_join(d_quantile_summary, boots_summary)

# plot summary plots of each habitat designation
pdf(paste('sequencing_rpoB/plots/analyses/habitat_preference/bootstrap_histograms_', otu_similarity, '.pdf', sep = ''), width = 8, height = 10)

# plot a couple of these
filter(boots2, habitat_preference == 'freshwater:mud_and_shore:terrestrial') %>%
  ggplot(., aes(prop_present)) +
  facet_grid(otu ~clust) +
  geom_histogram(col = 'black', fill = 'white', bins = 30) +
  geom_point(aes(quantile_threshold,10), shape = 21, fill = 'grey', size = 3) +
  theme_bw(base_size = 8) +
  geom_vline(aes(xintercept = prop_available), d_habitats, col = 'red') +
  labs(title = 'Generalists')
  

filter(boots2, habitat_preference == 'mud_and_shore:terrestrial') %>%
  ggplot(., aes(prop_present)) +
  facet_grid(otu ~clust) +
  geom_histogram(col = 'black', fill = 'white', bins = 30) +
  geom_point(aes(quantile_threshold,10), shape = 21, fill = 'grey', size = 3) +
  theme_bw(base_size = 8) +
  geom_vline(aes(xintercept = prop_available), d_habitats, col = 'red') +
  labs(title = 'Mud and Shore / Terrestrial specialists')

filter(boots2, habitat_preference == 'freshwater:mud_and_shore') %>%
  group_by(otu) %>%
  sample_n_groups(15) %>%
  ggplot(., aes(prop_present)) +
  facet_grid(otu ~clust) +
  geom_histogram(col = 'black', fill = 'white', bins = 30) +
  geom_point(aes(quantile_threshold,10), shape = 21, fill = 'grey', size = 3) +
  theme_bw(base_size = 8) +
  geom_vline(aes(xintercept = prop_available), d_habitats, col = 'red') +
  labs(title = 'Mud and Shore / Freshwater specialists')

filter(boots2, habitat_preference == 'freshwater:terrestrial') %>%
  group_by(otu) %>%
  sample_n_groups(15) %>%
  ggplot(., aes(prop_present)) +
  facet_grid(otu ~clust) +
  geom_histogram(col = 'black', fill = 'white', bins = 30) +
  geom_point(aes(quantile_threshold,10), shape = 21, fill = 'grey', size = 3) +
  theme_bw(base_size = 8) +
  geom_vline(aes(xintercept = prop_available), d_habitats, col = 'red') +
  labs(title = 'Freshwater / Terrestrial specialists')

filter(boots2, habitat_preference == 'freshwater') %>%
  group_by(otu) %>%
  sample_n_groups(15) %>%
  ggplot(., aes(prop_present)) +
  facet_grid(otu ~clust) +
  geom_histogram(col = 'black', fill = 'white', bins = 30) +
  geom_point(aes(quantile_threshold,10), shape = 21, fill = 'grey', size = 3) +
  theme_bw(base_size = 8) +
  geom_vline(aes(xintercept = prop_available), d_habitats, col = 'red') +
  labs(title = 'Freshwater specialists')

filter(boots2, habitat_preference == 'terrestrial') %>%
  group_by(otu) %>%
  sample_n_groups(15) %>%
  ggplot(., aes(prop_present)) +
  facet_grid(otu ~clust) +
  geom_histogram(col = 'black', fill = 'white', bins = 30) +
  geom_point(aes(quantile_threshold,10), shape = 21, fill = 'grey', size = 3) +
  theme_bw(base_size = 8) +
  geom_vline(aes(xintercept = prop_available), d_habitats, col = 'red') +
  labs(title = 'Terrestrial specialists')

filter(boots2, habitat_preference == 'mud_and_shore') %>%
  group_by(otu) %>%
  sample_n_groups(15) %>%
  ggplot(., aes(prop_present)) +
  facet_grid(otu ~clust) +
  geom_histogram(col = 'black', fill = 'white', bins = 30) +
  geom_point(aes(quantile_threshold,10), shape = 21, fill = 'grey', size = 3) +
  theme_bw(base_size = 8) +
  geom_vline(aes(xintercept = prop_available), d_habitats, col = 'red') +
  labs(title = 'Mud and Shore specialists')

dev.off()

# combine habitat preference datasets
d_pref_multiple2 <- d_pref_multiple %>%
  left_join(., boots_summary) %>%
  select(., otu, num_present, average_prop, average_abundance, habitat_preference, clust)

d_pref_all <- bind_rows(d_pref_single, d_pref_multiple2)

ggplot(d_pref_all, aes(num_present)) +
  facet_wrap(~clust) +
  geom_bar(aes(fill = habitat_preference), col = 'black') +
  theme_bw()

# does proportion broadly correlate with number present
ggplot(d_pref_all, aes(average_prop, average_abundance)) +
  geom_point(aes(col =habitat_preference)) +
  stat_smooth(method = 'lm') +
  facet_wrap(~clust) +
  theme_bw()

d_pref_all_summary <- group_by(d_pref_all, otu, habitat_preference) %>%
  summarise(average_prop = mean(average_prop),
            average_abundance = mean(average_abundance),
            habitats_present = n(),
            num_present = sum(num_present),
            .groups = 'drop') %>%
  mutate(habitat_preference_n = str_count(habitat_preference, ':') +1)

# does proportion broadly correlate with number present
ggplot(d_pref_all_summary, aes(num_present, log10(average_abundance))) +
  geom_point(aes(col=as.character(habitats_present)),alpha = 0.2) +
  stat_smooth(method = 'lm', se = FALSE, size = 0.9) +
  facet_wrap(~habitat_preference, ncol = 4) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.9, 0.2)) +
  palettetown::scale_colour_poke('Number of habitats\nOTU is present in', pokemon = 'lapras', spread = 3) +
  guides(color = guide_legend(override.aes = list(size = 6, alpha = 1))) +
  labs(x = 'Prevalence',
       y = 'Average log10 Abundance') +
  scale_x_continuous(n.breaks = 10)

#ggsave('sequencing_rpoB/plots/analyses/presence_vs_abundance.pdf', last_plot(), width = 12, height = 6)

group_by(d_pref_all_summary, habitat_preference, habitat_preference_n, habitats_present) %>%
  tally() %>%
  arrange(habitat_preference_n)

#------------------------#
# save everything out ####
#------------------------#

# save out number of otus
tibble(otu_similarity = otu_similarity, raw_otu_n = raw_taxa_n, prev_otu_n = prev_taxa_n) %>%
  write.csv(., paste('sequencing_rpoB/data/summary_stats/total_species_stats_', otu_similarity, '.csv', sep = ''), row.names = FALSE)

# save out habitat preferences
write.csv(d_pref_all_summary, paste('sequencing_rpoB/data/habitat_preference/habitat_preference_', otu_similarity, '.csv', sep = ''), row.names = FALSE)

# sample stats
# d_sample is first
d_ps2 <- filter(d_ps, abundance > 0) %>%
  select(sample, otu) %>%
  left_join(., select(d_pref_all_summary, otu, habitat_preference)) %>%
  group_by(sample, habitat_preference) %>%
  tally() %>%
  pivot_wider(names_from = 'habitat_preference', values_from = 'n') %>%
  mutate(across(everything(), function(x) replace_na(x, 0))) %>%
  select(., sample, freshwater, terrestrial, mud_and_shore, `freshwater:terrestrial`, `freshwater:mud_and_shore`, `mud_and_shore:terrestrial`, `freshwater:mud_and_shore:terrestrial`)

left_join(d_sample, d_ps2) %>%
  write.csv(., paste('sequencing_rpoB/data/summary_stats/sample_stats_', otu_similarity, '.csv', sep = ''), row.names = FALSE)

# lets look at a generalist
d_habitats <- mutate(d_habitats, clust = ifelse(clust == 'mud_and_shore', 'marine mud', clust))



# plot different distributions for the presentation.
filter(boots2, habitat_preference == 'freshwater:terrestrial') %>%
  filter(otu == sample(otu, 1)) %>%
  mutate(., clust = ifelse(clust == 'mud_and_shore', 'marine mud', clust)) %>%
  ggplot(., aes(prop_present)) +
  facet_wrap(~clust) +
  geom_histogram(aes(fill = clust), col = 'black') +
  geom_point(aes(quantile_threshold,10), shape = 21, fill = 'grey', size = 8) +
  theme_bw(base_size = 24) +
  geom_vline(aes(xintercept = prop_available), d_habitats, col = 'red', size = 1.5) +
  theme(legend.position = "none",
        axis.text.y = element_blank()) +
  scale_fill_manual(values = c('#1E90FE', 'sienna4', 'green4')) +
  labs(x = 'Proportion used',
       y = '') +
  xlim(c(-0.1,1))

ggsave('sequencing_rpoB/plots/analyses/hab_pref_pres_4.pdf', last_plot(), width = 12, height = 4)




