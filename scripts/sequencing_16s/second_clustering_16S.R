#-------------------------------------------#
# look at clustering of otus and samples ####
#-------------------------------------------#

# decisions we made
# remove estuarine mud low polyhaline, thrift rhizosphere, and beach supratidal from clustering analysis

# what does this script do

# 1. runs PCoA on filtered sample set
# 2. runs medoid clustering approach on distance matrix
# 3. runs hierarchical clustering on distance matrix

# clean workspace
rm(list = ls())

# set seed to always get the same clustering assignments
set.seed(42)

# load packages ####
library(microViz) # remotes::install_github("david-barnett/microViz")
library(phyloseq)
library(vegan)
library(patchwork)
library(tidyverse)
library(MicrobioUoE)
library(cluster)
library(NbClust)
library(fpc)
library(dendextend)
library(factoextra)
remotes::install_github("leffj/mctoolsr")

path_fig <- 'sequencing_16S/plots'
path_fig <- 'plots/sequencing_16s'

# source extra functions
source('scripts/sequencing_16S/extra_functions.R')

# load data
ps <- readRDS('data/sequencing_16s/ps_16s_low_depth_removed.rds')

meta <- sample_data(ps) %>% data.frame()

#-----------------------------------------------------------------------#
# 1. Remove habitats with low sample size/too similar and re-do PCoA ####
#-----------------------------------------------------------------------#

# remove some habitats we do not want
to_keep <- filter(meta, !habitat_group_16s %in% c('estuarine mud_low polyhaline', 'thrift_rhizosphere', 'beach_supratidal')) %>%
  row.names(.)

ps <- prune_samples(to_keep, ps)
meta <- sample_data(ps) %>% data.frame()
summary <- group_by(meta, habitat_group_16s) %>%
  tally()

# transform counts to relative abundances for ordination
ps_prop <- transform_sample_counts(ps, function(x){x / sum(x)})

# wrangle the metadata
d_samp <- data.frame(sample_data(ps_prop))
d_samp <- mutate(d_samp, location_fac = as.factor(location),
                 habitat_group_fac = as.factor(habitat_group_16s)) %>%
  unite(., 'id', location_fac, habitat_group_fac, sep = ':', remove = FALSE)

# calculate distance matrix
ps_wunifrac <- phyloseq::distance(ps_prop, method = 'wunifrac')

# run a betadisper
mod_betadisper <- betadisper(ps_wunifrac, d_samp$habitat_group_fac)

# grab centroids and other data
d_fig <- get_betadisper_data(mod_betadisper)

# combine centroid and eigenvector dataframes for plotting
betadisper_lines <- merge(select(d_fig$centroids, group, PCoA1, PCoA2), select(d_fig$eigenvector, group, PCoA1, PCoA2, sample), by = c('group'))

# add distances to eigenvector and lines data
betadisper_lines <- mutate(betadisper_lines, distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))
d_fig$eigenvector$distances <- d_fig$distances$distances

# split up group into week and plastic context
betadisper_lines <- left_join(betadisper_lines, rownames_to_column(meta, var = 'sample') %>% select(sample, location))
d_fig$eigenvector <- left_join(d_fig$eigenvector, rownames_to_column(meta, var = 'sample') %>% select(sample, location)) 

d_fig$eigenvalue <- mutate(d_fig$eigenvalue, percent = eig/sum(eig)*100)

# get correct eigenvalue
correct_eigenvalues <- ape::pcoa(ps_wunifrac, correction = 'cailliez', d_samp$id) %>%
  .$values %>%
  pull(Rel_corr_eig)
correct_eigenvalues[1:4]

# plot scree plot
ggplot(tibble(eig = correct_eigenvalues, n = 1:length(correct_eigenvalues)), aes(n, eig)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  labs(y = 'relative corrected eigenvalue',
       x = 'PCoA axis')

# set up colours based on habitats
cols <- tibble(group = c("woodland_oak", "estuarine mud_low polyhaline", "woodland_pine", "reservoir", "river", "estuarine mud_oligohaline", "estuarine mud_full saline", "beach_supratidal", "pasture", "beach_subtidal","thrift_rhizosphere","beach_seaweed","field_wheat","rock_samphire","marine mud_full saline"),
               col = c('#089a2d', '#995a08', '#106c12', '#1170bd', '#9dcdf4', '#663c05', '#b6966b', '#f5e279', '#61dd1e', '#f2f426', '#c5f8ae', '#a6ab52', '#9ff121', '#5e8128', '#714a03'),
               hab_order = c(1.1, 2.1, 1.2, 3.1, 3.2, 2.2, 2.3, 2.4, 1.3, 2.5, 2.6, 2.7, 1.4, 1.5, 2.8))
cols <- filter(cols, group %in% d_samp$habitat_group_16s)
cols <- mutate(cols, habitat_group_16s = group) %>% arrange(hab_order)

# change colours names
cols <- mutate(cols, group2 = gsub('_', ' ', group))

# plot the samples across the first 4 axes
d_fig$eigenvector %>%
  select(group:PCoA4) %>%
  mutate(to_order = PCoA1) %>%
  pivot_longer(contains('PCoA'), names_to = 'axis', names_prefix = 'PCoA', values_to = 'eigenvector') %>%
  mutate(habitat2 = gsub('_', '\n', group),
         axis = case_when(axis == 1 ~ 'Axis 1: 34.69%',
                          axis == 2 ~ 'Axis 2: 14.24%',
                          axis == 3 ~ 'Axis 3: 8.32%',
                          axis == 4 ~ 'Axis 4: 5.96%')) %>%
  ggplot(., aes(forcats::fct_reorder(habitat2, to_order), eigenvector)) +
  geom_hline(aes(yintercept = 0)) +
  geom_point(aes(col = group), size = 2) +
  facet_wrap(~axis) +
  scale_color_manual('Habitat', values = setNames(cols$col, cols$group2)) +
  labs(x = 'Habitat',
       y = 'Eigenvector') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.7)) 

ggsave(file.path(path_fig, 'ordination_axis_split.png'), last_plot(), width = 12, height = 8)
ggsave(file.path(path_fig, 'ordination_axis_split.pdf'), last_plot(), width = 12, height = 8)

# plot PCoA
p1 <- ggplot() +
  geom_segment(aes(x = PCoA1.x*-1, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y*-1, group = sample, col = group), betadisper_lines, show.legend = FALSE) +
  geom_point(aes(PCoA1*-1, PCoA2, col = group, shape = location),d_fig$eigenvector) +
  geom_point(aes(PCoA1*-1, PCoA2, fill = group), d_fig$centroids, size = 5, col = 'white', shape = 21, show.legend = FALSE) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12),
        legend.position = 'right') +
  scale_color_manual('Habitat', values = setNames(cols$col, cols$group), labels = function(x) gsub("_", " ", x)) +
  scale_fill_manual('Habitat', values = setNames(cols$col, cols$group), labels = function(x) gsub("_", " ", x)) +
  scale_shape_discrete('Location') +
  labs(y = 'Axis 2 (14.24%)',
       x = 'Axis 1 (34.69%)')

# save out plot
saveRDS(p1, file.path(path_fig, 'ordination_16S_clean.rds'))

ggsave(file.path(path_fig, 'ordination_16S_clean.pdf'), p1, width = 10, height = 7)
ggsave(file.path(path_fig, 'ordination_16S_clean.png'), p1, width = 10, height = 7)

# look at pairwise differences between groups
# pairwise permanovas
mods <- calc_pairwise_permanovas(ps_wunifrac, d_samp, 'habitat_group_fac', n_perm = 9999)
filter(mods, pvalFDR >= 0.05) %>%
  select(., X1, X2, R2, pval, pvalFDR)
# all significantly different now

#---------------------------------------------------#
# 2. cluster 16S samples using medoid clustering ####
#---------------------------------------------------#

# this follows the tutorials in Modern Statistics for Modern Biology and references therein
# https://www.huber.embl.de/msmb/Chap-Clustering.html

# use pam clustering which uses the uses medoid centreing
# use the package cluster::pam
# run clusGap to look at which number of clusters is best
# also run NbClust

# choose the number of clusters
pamfun = function(x, k){list(cluster = pam(x, k, cluster.only = TRUE))}

# got to use a matrix here
# use the scores for the PCoA of the plots - taken from phyloseq vignette
d_pcoa_samples <- dplyr::select(d_fig$eigenvector, sample, starts_with('PCoA')) %>%
  column_to_rownames('sample')
pam_clusters = clusGap(d_pcoa_samples, FUN = pamfun, K.max = 12, B = 300, verbose = TRUE)
d_gap <- data.frame(pam_clusters$Tab, k=1:nrow(pam_clusters$Tab)) %>%
  data.frame()

# this errors!
#pam_nbclust <- NbClust(d_pcoa_samples, method = 'kmeans', max.nc = 12)
# because of negative eigenvectors

# create a pcoa data set - choose only eigenvalues that are > 0
#https://stackoverflow.com/questions/8924488/applying-the-pvclust-r-function-to-a-precomputed-dist-object#27148408
# error message tells how may eigenvalues are > 0
d_pcoa_correct <- cmdscale(ps_wunifrac, length(labels(ps_wunifrac))-1, eig=TRUE)
d_pcoa_correct <- cmdscale(ps_wunifrac, sum(d_pcoa_correct$eig > 0))

# redo clustering on corrected PCoA plot
pam_clusters = clusGap(d_pcoa_correct, FUN = pamfun, K.max = 12, B = 300, verbose = TRUE)
d_gap <- data.frame(pam_clusters$Tab, k=1:nrow(pam_clusters$Tab)) %>%
  data.frame()

# calculate the optimal number of clusters
maxSE(d_gap$gap, d_gap$SE.sim, method = 'Tibs2001SEmax')

# do NbClust on the pcoa_correct
pam_nbclust <- NbClust(d_pcoa_correct, method = 'kmeans', max.nc = 12) # uses kmeans instead of k-medoids
# says the best is 3

# visualise both cluster sets
clustering <- pam(d_pcoa_correct, k = 3)

# put cluster assignment into dataframe
# align cluster assignments when numbers do not match
cluster_numbers <- data.frame(sample = names(clustering$clustering), cluster_gap = clustering$clustering,
                              cluster_nbclust = pam_nbclust$Best.partition) %>%
  mutate(across(starts_with('cluster'), as.character),
         cluster_nbclust = case_when(cluster_nbclust == '2' ~ '1',
                                     cluster_nbclust == '3' ~ '2',
                                     cluster_nbclust == '1' ~ '3')) %>%
  pivot_longer(starts_with('cluster'), names_to = 'method', names_prefix = 'cluster_', values_to = 'cluster')

d_clustering <- rownames_to_column(d_pcoa_samples, var = 'sample') %>%
  left_join(cluster_numbers)

# make dataset for medoid clustering
cluster_medoid <- data.frame(sample = names(clustering$clustering), medoid_gap = clustering$clustering,
                             medoid_nbclust = pam_nbclust$Best.partition) %>%
  mutate(across(starts_with('cluster'), as.character),
         medoid_nbclust = case_when(medoid_nbclust == '2' ~ '1',
                                    medoid_nbclust == '3' ~ '2',
                                    medoid_nbclust == '1' ~ '3'))

d_centroids <- d_clustering %>%
  group_by(., cluster, method) %>%
  summarise(across(starts_with('PCoA'), mean), .groups = 'drop')

d_lines <- merge(select(d_clustering, sample, method, cluster, PCoA1, PCoA2), select(d_centroids, method, cluster, PCoA1, PCoA2), by = c('method', 'cluster')) %>%
  mutate(distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))

p3 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = cluster),d_clustering, show.legend = TRUE) +
  geom_point(aes(PCoA1, PCoA2, col = cluster), d_centroids, size = 3, show.legend = FALSE) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = sample, col = cluster), d_lines, show.legend = FALSE) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12)) +
  facet_wrap(~method)

ggsave(file.path(path_fig, 'partition_clustering.png'), p3, width = 10, height = 5)

# plot the gap statistic
ggplot(d_gap, aes(k, gap)) +
  geom_line() +
  geom_point(size = 5) +
  geom_linerange(aes(ymin = gap - SE.sim, ymax = gap + SE.sim)) +
  theme_bw(base_size = 14) +
  labs(y = 'Gap statistic',
       x = 'Number of clusters') +
  scale_x_continuous(breaks = 1:12)

rownames(d_samp) == rownames(cluster_numbers)

#---------------------------------------------------------#
# 3. cluster 16S samples using hierarchical clustering ####
#---------------------------------------------------------#

# do clustering using different methods
hier_clust <- agnes(ps_wunifrac, method = 'ward')
hier_clust2 <- agnes(ps_wunifrac, method = 'average')
hier_clust3 <- agnes(ps_wunifrac, method = 'complete')
hier_clust4 <- agnes(ps_wunifrac, method = 'single')

# look at the agglomerative coefficient
hier_clust$ac
hier_clust2$ac
hier_clust3$ac
hier_clust4$ac
# ward method gives the strongest clusters based on agglomerative coefficienct

# look at the correlation between clustering and the actual data
cor(cophenetic(hier_clust), ps_wunifrac)
cor(cophenetic(hier_clust2), ps_wunifrac)
cor(cophenetic(hier_clust3), ps_wunifrac)
cor(cophenetic(hier_clust4), ps_wunifrac)
# average method gives the best representation of the data

# calculate the optimum number of clusters as earlier
# do this on both average and ward methods and compare
hier_clusgap <- clusGap(d_pcoa_correct, FUN = factoextra::hcut, K.max = 12, B = 300, hc_func = 'agnes', hc_method = 'ward')
hier_clusgap2 <- clusGap(d_pcoa_correct, FUN = factoextra::hcut, K.max = 12, B = 300, hc_func = 'agnes', hc_method = 'average')

d_gap <- data.frame(hier_clusgap$Tab, k=1:nrow(hier_clusgap$Tab)) %>%
  data.frame()
d_gap2 <- data.frame(hier_clusgap2$Tab, k=1:nrow(hier_clusgap2$Tab)) %>%
  data.frame()

# calculate the optimal number of clusters
maxSE(d_gap$gap, d_gap$SE.sim, method = 'firstSEmax')
maxSE(d_gap$gap, d_gap$SE.sim, method = 'Tibs2001SEmax')
maxSE(d_gap2$gap, d_gap2$SE.sim, method = 'firstSEmax')
maxSE(d_gap2$gap, d_gap2$SE.sim, method = 'Tibs2001SEmax')
# average says there is only a single cluster - not right

# decided to use the ward method, but hard to test which is "better"

# can compare the dendograms of different clusterings

# compare dendograms using dendextend
dend_ward <- as.dendrogram(hier_clust)
dend_average <-as.dendrogram(hier_clust2)
dend_single <- as.dendrogram(hier_clust4)
dend_complete <- as.dendrogram(hier_clust3)
list_of_dends <- dendlist(dend_ward, dend_average, dend_single, dend_complete)
names(list_of_dends) <- c('ward', 'average', 'single', 'complete')

# calculate correlation between dendograms
cor.dendlist(list_of_dends)
cor.dendlist(list_of_dends, method = "common")

# ward, average and complete are all pretty similar

# calculate entanglement between dendograms
# 0 = no entanglement, 1 = full entanglement
entanglement(dend_ward, dend_average)
entanglement(dend_ward, dend_complete)
entanglement(dend_ward, dend_single)

# do NbClustering on Ward method
hier_nbclust <- NbClust::NbClust(d_pcoa_correct, method = 'ward.D', max.nc = 12)
summary(hier_nbclust)

ggplot(d_gap2, aes(k, gap)) +
  geom_line() +
  geom_point(size = 5) +
  geom_linerange(aes(ymin = gap - SE.sim, ymax = gap + SE.sim)) +
  theme_bw(base_size = 14) +
  labs(y = 'Gap statistic',
       x = 'Number of clusters') +
  scale_x_continuous(breaks = 1:12)

# assign clusters for hierarchical clusterings 
clustering <- cutree(hier_clust, k = 3)
cluster_numbers <- data.frame(sample = names(clustering), cluster_gap = clustering,
                              cluster_nbclust = hier_nbclust$Best.partition) %>%
  mutate(across(starts_with('cluster'), as.character)) %>%
  pivot_longer(starts_with('cluster'), names_to = 'method', names_prefix = 'cluster_', values_to = 'cluster')

cluster_hclust <- data.frame(sample = names(clustering), hier_gap = clustering,
                             hier_nbclust = hier_nbclust$Best.partition) %>%
  mutate(across(starts_with('cluster'), as.character))

d_clustering <- rownames_to_column(d_pcoa_samples, var = 'sample') %>%
  left_join(cluster_numbers)
d_centroids <- d_clustering %>%
  group_by(., cluster, method) %>%
  summarise(across(starts_with('PCoA'), mean), .groups = 'drop')
d_lines <- merge(select(d_clustering, sample, method, cluster, PCoA1, PCoA2), select(d_centroids, method, cluster, PCoA1, PCoA2), by = c('method', 'cluster')) %>%
  mutate(distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))

p5 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = cluster),d_clustering, show.legend = FALSE) +
  geom_point(aes(PCoA1, PCoA2, col = cluster), d_centroids, size = 3, show.legend = FALSE) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = sample, col = cluster), d_lines, show.legend = FALSE) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12)) +
  facet_wrap(~method)

ggsave(file.path(path_fig, 'hierarchical_clustering_groups.png'), p5, width = 10, height = 5)

# create dataframe for cluster assignments for each sample
d_clusters <- left_join(cluster_medoid, cluster_hclust) %>%
  mutate(across(everything(), as.character))

write.csv(d_clusters, 'data/sequencing_16S/sample_cluster_assignments.csv', row.names = FALSE)

# compare clusters
d_compare_cluster <- left_join(d_clusters, rownames_to_column(d_samp, var = 'sample')) %>%
  pivot_longer(cols = c(medoid_gap:hier_nbclust), names_to = 'method', values_to = 'cluster') %>%
  group_by(habitat_group_16s, cluster, method) %>%
  tally() %>%
  ungroup() %>%
  mutate(habitat2 = gsub('_', '\n', habitat_group_16s)) 

ggplot(d_compare_cluster, aes(x = forcats::fct_reorder(habitat2, cluster), y = n, fill = cluster)) +
  geom_bar(stat = 'identity', col = 'black', show.legend = FALSE) +
  theme_bw(base_size = 14) +
  labs(y = 'Number of samples',
       x = 'habitat type') +
  facet_wrap(~method) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(path_fig, 'partition_cluster_assignment.png'), last_plot(), width = 11, height = 9)
# partition assignment here is great

# make final clustering plot 

# plot out the medoid nbclust samples because they are the ones we observed

# make the clustered datasets
d_clustering <- rownames_to_column(d_pcoa_samples, var = 'sample') %>%
  left_join(cluster_medoid)
d_centroids <- d_clustering %>%
  group_by(., medoid_nbclust) %>%
  summarise(across(starts_with('PCoA'), mean), .groups = 'drop')

d_lines <- merge(select(d_clustering, sample, medoid_nbclust, PCoA1, PCoA2), select(d_centroids, medoid_nbclust, PCoA1, PCoA2), by = c('medoid_nbclust')) %>%
  mutate(distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))

# make the plot

ggplot() +
  geom_segment(aes(x = PCoA1.x*-1, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y*-1, group = sample, col = medoid_nbclust), d_lines, show.legend = FALSE) +
  geom_point(aes(PCoA1*-1, PCoA2, fill = medoid_nbclust), shape = 21, col = 'white', d_clustering) +
  geom_point(aes(PCoA1*-1, PCoA2, fill = medoid_nbclust), shape = 21, col = 'white', d_centroids, size = 5, show.legend = FALSE) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12),
        legend.position = 'right') +
  scale_fill_manual('Biome', values = c('#53C20A', '#5ECAE2', '#3911EE'), labels = c('Land', 'Freshwater', 'Marine')) +
  scale_color_manual('Habitat cluster', values = c('#53C20A', '#5ECAE2', '#3911EE'), labels = c('Land', 'Freshwater', 'Marine')) +
  labs(y = 'Axis 2 (14.24%)',
       x = 'Axis 1 (34.69%)')

# save out object
saveRDS(last_plot(), file.path(path_fig, 'cluster_pcoa.rds'))