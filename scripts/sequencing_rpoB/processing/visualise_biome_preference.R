#--------------------------#
# what this script does ####
#--------------------------#

# visualises some of the biome preferences as assigned by the bootstrapping of presence/absences

# load in necessary packages
librarian::shelf(phyloseq, tidyverse)

# read in the asv habitat preference data
hab_pref_boots <- readRDS('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/bootstraps/habpref_boots_asv.rds')

hab_pref <- read.csv('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_asv.csv')

# read in clusters
# read in cluster assignments
clusters <- read.csv('data/sequencing_16S/sample_cluster_assignments.csv') %>%
  mutate(clust = case_when(medoid_nbclust == '2' ~ 'freshwater',
                           medoid_nbclust == '3' ~ 'marine',
                           medoid_nbclust == '1' ~ 'land'),
         clust = ifelse(sample == 'sample_s46', 'marine', clust),
         clust_fac = as.factor(clust)) %>%
  select(., sample, clust)

# look at availability of clusters
d_habitats <- clusters %>% 
  select(sample, clust) %>%
  distinct() %>%
  group_by(clust) %>%
  tally() %>%
  mutate(prop_available = n / sum(n)) %>%
  rename(num_available = n)

# look at quantiles of habitat preference
# calculate which quantile the available proportion is in of the distribution of used proportions
d_quantiles <- group_by(hab_pref_boots, otu, clust) %>%
  summarise(quantile = sum(prop_present <= unique(prop_available))/n(), .groups = 'drop',
            quantile_threshold = quantile(prop_present, 0.975),
            prop_available = unique(prop_available))

# make a plot of distributions of one of each of the habitat preferences

to_plot <- group_by(hab_pref, habitat_preference) %>%
  filter(habitats_present > 1) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(habitat_preference = case_when(habitat_preference == 'freshwater' ~ 'freshwater\nspecialist',
                                        habitat_preference == 'marine_mud' ~ 'marine\nspecialist',
                                        habitat_preference == 'terrestrial' ~ 'land\nspecialist',
                                        habitat_preference == 'freshwater:marine_mud' ~ 'freshwater + marine\ngeneralist',
                                        habitat_preference == 'freshwater:terrestrial' ~ 'freshwater + land\ngeneralist',
                                        habitat_preference == 'marine_mud:terrestrial' ~ 'marine + land\ngeneralist',
                                        habitat_preference == 'freshwater:marine_mud:terrestrial' ~ 'full\ngeneralist'))

hab_pref_boots %>%
  filter(otu %in% to_plot$otu) %>%
  left_join(., select(to_plot, otu, habitat_preference)) %>%
  mutate(clust = case_when(clust == 'terrestrial' ~ 'land',
                           clust == 'marine_mud' ~ 'marine',
                           TRUE ~ 'freshwater')) %>%
  ggplot(., aes(prop_present)) +
  facet_grid(habitat_preference ~ clust) +
  geom_histogram(aes(fill = clust), col = 'black', bins = 40, show.legend = FALSE) +
  theme_bw(base_size = 10) +
  geom_vline(aes(xintercept = prop_available), d_habitats, col = 'black') +
  labs(y = 'Number of 1000 bootstraps',
       x = 'Proportional biome presence') +
  scale_fill_manual('Biome', values = c('#5ECAE2', '#53C20A', '#3911EE'))

ggsave('plots/manuscript_plots/visualise_habitat_preference.png', last_plot(), width = 8, height = 8.5)
