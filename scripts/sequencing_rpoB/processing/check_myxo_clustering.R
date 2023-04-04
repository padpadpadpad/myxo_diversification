#-----------------------------------------------#
# check the clustering of the myxo rpoB data ####
#-----------------------------------------------#

# what does this script do

# 1. run PCoA on all prevalence filtered myxo samples based on habitat
# 2. run pairwise permanovas on all habitats to see which are not significantly different

# clean workspace
rm(list = ls())

# load packages ####
librarian::shelf(david-barnett/microViz, phyloseq, vegan, patchwork, tidyverse, padpadpadpad/MicrobioUoE, ape)

# source extra functions
source('scripts/sequencing_16s/extra_functions.R')

# load data
ps <- readRDS('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_asv_filt.rds')
ps

# load in tree
tree <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_chronopl10.tre')

# add tree in
phy_tree(ps) <- tree

ps

# transform counts to relative abundances for ordination
ps_prop <- transform_sample_counts(ps, function(x){x / sum(x)})

# read in cluster assignments
clusters <- read.csv('data/sequencing_16S/sample_cluster_assignments.csv')

#-----------------------------------------------------------------------#
# 1. Run PCoA looking at community composition of different habitats ####
#-----------------------------------------------------------------------#

# wrangle the metadata
d_samp <- data.frame(sample_data(ps_prop)) %>%
  rownames_to_column(var = 'sample') %>%
  left_join(., clusters)
d_samp <- mutate(d_samp, location_fac = as.factor(location),
                 habitat_group_fac = as.factor(habitat_group_16s)) %>%
  mutate(clust = case_when(medoid_nbclust == '2' ~ 'freshwater',
                           medoid_nbclust == '3' ~ 'mud_and_shore',
                           medoid_nbclust == '1' ~ 'terrestrial'),
         clust_fac = as.factor(clust))

# calculate distance matrix
ps_wunifrac <- phyloseq::distance(ps_prop, method = 'wunifrac')

# run a betadisper
mod_betadisper <- betadisper(ps_wunifrac, d_samp$clust_fac)

plot(mod_betadisper)

# grab centroids and other data
d_fig <- get_betadisper_data(mod_betadisper)

# combine centroid and eigenvector dataframes for plotting
betadisper_lines <- merge(select(d_fig$centroids, group, PCoA1, PCoA2), select(d_fig$eigenvector, group, PCoA1, PCoA2, sample), by = c('group'))

# add distances to eigenvector and lines data
betadisper_lines <- mutate(betadisper_lines, distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))
d_fig$eigenvector$distances <- d_fig$distances$distances

# split up group into week and plastic context
betadisper_lines <- left_join(betadisper_lines, select(d_samp, sample, clust))
d_fig$eigenvector <- left_join(d_fig$eigenvector, select(d_samp, sample, clust)) 

d_fig$eigenvalue <- mutate(d_fig$eigenvalue, percent = eig/sum(eig)*100)

# get correct eigenvalue
correct_eigenvalues <- ape::pcoa(ps_wunifrac, correction = 'cailliez', d_samp$id) %>%
  .$values %>%
  pull(Rel_corr_eig)
correct_eigenvalues[1:2]

# plot PCoA
p1 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = group),d_fig$eigenvector) +
  geom_point(aes(PCoA1, PCoA2, col = group), d_fig$centroids, size = 5, show.legend = FALSE) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = sample, col = group), betadisper_lines) +
  #theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12),
        legend.position = 'right') +
  labs(y = 'Axis 2 (10.58%)',
       x = 'Axis 1 (28.22%)') +
  ggrepel::geom_label_repel(aes(PCoA1, PCoA2, col = group, label = sample),d_fig$eigenvector)

p1
