# look at clustering of otus and samples ####

# clean workspace
rm(list = ls())

# load packages ####
library(microViz)
library(DESeq2) # BiocManager::install("DESeq2")
library(phyloseq)
library(vegan)
library(patchwork) 
library(lme4)
library(flextable)
library(tidyverse)
library(palettetown)
library(MicrobioUoE)
library(clusterMany)
library(cluster)
library(ggtree)
# library(mclust)
library(NbClust)
library(fpc)
library(dendextend)
library(factoextra)
library(scclust)

# if not installed, install mctoolsr run remotes::install_github('leffj/mctoolsr')

path_fig <- 'sequencing_16S/plots/analyses'

# source extra functions
source('sequencing_16S/scripts/extra_functions.R')

# load data
#ps <- readRDS('data/sequencing/output/run_4/ps_rarefied.rds')
ps <- readRDS('sequencing_16S/data/output/run_merged_runs_new/ps_low_depth_removed_otu97.rds')

# replace with new metadata
meta <- read.csv('sequencing_16S/data/metadata_complete_fixed.csv', stringsAsFactors = FALSE)
meta <- mutate(meta, habitat_group2 = ifelse(habitat_group %in% c('beach_supratidal', 'thrift_rhizosphere'), 'beach_thrift', habitat_group))
row.names(meta) <- paste('sample_s', meta$id, sep = '')
sample_data(ps) <- sample_data(meta)

# remove some habitats we do not want
to_keep <- filter(meta, !habitat_group %in% c('estuarine mud_low polyhaline', 'thrift_rhizosphere', 'beach_supratidal')) %>%
  row.names(.)

ps <- prune_samples(to_keep, ps)

meta <- sample_data(ps) %>% data.frame()

summary <- group_by(meta, habitat_group2) %>%
  tally()

sort(sample_sums(ps))
rank_names(ps)

# transform counts to relative abundances for ordination
ps_prop <- transform_sample_counts(ps, function(x){x / sum(x)})

# wrangle the metadata
d_samp <- data.frame(sample_data(ps_prop))
d_samp <- mutate(d_samp, location_fac = as.factor(location),
                 habitat_group_fac = as.factor(habitat_group)) %>%
  unite(., 'id', location_fac, habitat_group_fac, sep = ':', remove = FALSE)

# calculate distance matrix
ps_wunifrac <- phyloseq::distance(ps_prop, method = 'wunifrac')

# to label - things that are on average far away from each other
to_label <- dist_2_df(ps_wunifrac) %>%
  pivot_longer(., starts_with('X'), names_to = 'id', values_to = 'sample_id') %>%
  group_by(sample_id) %>%
  summarise(dist = median(dist),
            .groups = 'drop') %>%
  filter(dist > quantile(dist, 0.95))

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

ggplot(tibble(eig = correct_eigenvalues, n = 1:length(correct_eigenvalues)), aes(n, eig)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  labs(y = 'relative corrected eigenvalue',
       x = 'PCoA axis')

# set up colours based on habitats
cols <- tibble(group = c("woodland_oak", "estuarine mud_low polyhaline", "woodland_pine", "reservoir", "river", "estuarine mud_oligohaline", "estuarine mud_full saline", "beach_thrift", "pasture", "beach_subtidal","thrift_rhizosphere","beach_seaweed","field_wheat","rock_samphire","marine mud_full saline"),
               col = c('#089a2d', '#995a08', '#106c12', '#1170bd', '#9dcdf4', '#663c05', '#b6966b', '#f5e279', '#61dd1e', '#f2f426', '#c5f8ae', '#a6ab52', '#9ff121', '#5e8128', '#714a03'))
cols <- filter(cols, group %in% d_samp$habitat_group)
cols <- mutate(cols, habitat_group = group)

# plot the samples across the first 4 axes
d_fig$eigenvector %>%
  select(group:PCoA4) %>%
  mutate(to_order = PCoA1) %>%
  pivot_longer(contains('PCoA'), names_to = 'axis', names_prefix = 'PCoA', values_to = 'eigenvector') %>%
  mutate(habitat2 = gsub('_', '\n', group),
         axis = case_when(axis == 1 ~ 'Axis 1: 39.97%',
                          axis == 2 ~ 'Axis 2: 14.01%',
                          axis == 3 ~ 'Axis 3: 7.72%',
                          axis == 4 ~ 'Axis 4: 5.12%')) %>%
  ggplot(., aes(forcats::fct_reorder(habitat2, to_order), eigenvector)) +
  geom_hline(aes(yintercept = 0)) +
  geom_point(aes(col = group), size = 2) +
  facet_wrap(~axis) +
  scale_color_manual('Habitat', values = setNames(cols$col, cols$group)) +
  labs(x = 'Habitat',
       y = 'Eigenvector') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(file.path(path_fig, 'ordination_axis_split.png'), last_plot(), width = 12, height = 8)
ggsave(file.path(path_fig, 'ordination_axis_split.pdf'), last_plot(), width = 12, height = 8)

# plot PCoA
p1 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = group, shape = location),d_fig$eigenvector) +
  #ggrepel::geom_label_repel(aes(PCoA1, PCoA2, col = habitat, label = id), data = to_label) +
  geom_point(aes(PCoA1, PCoA2, col = group), d_fig$centroids, size = 3, show.legend = FALSE) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = sample, col = group), betadisper_lines) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12),
        legend.position = 'right') +
  scale_color_manual('Habitat', values = setNames(cols$col, cols$group)) +
  labs(x = 'Habitat',
       y = 'Eigenvector')

# look at pairwise differences between groups
# pairwise permanovas
mods <- calc_pairwise_permanovas(ps_wunifrac, d_samp, 'habitat_group_fac', n_perm = 9999)
filter(mods, pvalFDR >= 0.05) %>%
  select(., X1, X2, R2, pval, pvalFDR)

#--------------------------------------------#
# do clustering to try and group habitats ####
#--------------------------------------------#

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
pam_nbclust <- NbClust(d_pcoa_samples, method = 'kmeans', max.nc = 12)

# create a pcoa data set - choose only eigenvalues that are > 0
#https://stackoverflow.com/questions/8924488/applying-the-pvclust-r-function-to-a-precomputed-dist-object#27148408
d_pcoa_correct <- cmdscale(ps_wunifrac, 43)
pam_clusters = clusGap(d_pcoa_correct, FUN = pamfun, K.max = 12, B = 300, verbose = TRUE)
d_gap <- data.frame(pam_clusters$Tab, k=1:nrow(pam_clusters$Tab)) %>%
  data.frame()

# calculate the optimal number of clusters
maxSE(d_gap$gap, d_gap$SE.sim, method = 'firstSEmax')
maxSE(d_gap$gap, d_gap$SE.sim, method = 'Tibs2001SEmax')

pam_nbclust <- NbClust(d_pcoa_correct, method = 'kmeans', max.nc = 12) # uses kmeans instead of k-medoids
# says the best is 3

# see the cluster
clustering <- pam(d_pcoa_correct, k = 12)
cluster_numbers <- data.frame(sample = names(clustering$clustering), cluster_gap = clustering$clustering,
                              cluster_nbclust = pam_nbclust$Best.partition) %>%
  mutate(across(starts_with('cluster'), as.character)) %>%
  pivot_longer(starts_with('cluster'), names_to = 'method', names_prefix = 'cluster_', values_to = 'cluster')
cluster_medoid <- data.frame(sample = names(clustering$clustering), medoid_gap = clustering$clustering,
                             medoid_nbclust = pam_nbclust$Best.partition) %>%
  mutate(across(starts_with('cluster'), as.character))

d_clustering <- rownames_to_column(d_pcoa_samples, var = 'sample') %>%
  left_join(cluster_numbers)
d_centroids <- d_clustering %>%
  group_by(., cluster, method) %>%
  summarise(across(starts_with('PCoA'), mean), .groups = 'drop')
d_lines <- merge(select(d_clustering, sample, method, cluster, PCoA1, PCoA2), select(d_centroids, method, cluster, PCoA1, PCoA2), by = c('method', 'cluster')) %>%
  mutate(distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))

p3 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = cluster),d_clustering, show.legend = FALSE) +
  #ggforce::geom_mark_hull(aes(PCoA1, PCoA2, col = cluster), d_clustering, show.legend = FALSE, concavity = 2) +
  #ggrepel::geom_label_repel(aes(PCoA1, PCoA2, col = habitat, label = id), data = to_label) +
  geom_point(aes(PCoA1, PCoA2, col = cluster), d_centroids, size = 3, show.legend = FALSE) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = sample, col = cluster), d_lines, show.legend = FALSE) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12)) +
  facet_wrap(~method)

ggplot(d_gap, aes(k, gap)) +
  geom_line() +
  geom_point(size = 5) +
  geom_linerange(aes(ymin = gap - SE.sim, ymax = gap + SE.sim)) +
  theme_bw(base_size = 14) +
  labs(y = 'Gap statistic',
       x = 'Number of clusters') +
  scale_x_continuous(breaks = 1:12)

rownames(d_samp) == rownames(cluster_numbers)

d_compare_cluster <- left_join(cluster_numbers, rownames_to_column(d_samp, var = 'sample')) %>%
  group_by(habitat_group, cluster, method) %>%
  tally() %>%
  ungroup() %>%
  mutate(habitat2 = gsub('_', '\n', habitat_group))

ggplot(d_compare_cluster, aes(x = forcats::fct_reorder(habitat2, cluster), y = n, fill = cluster)) +
  geom_bar(stat = 'identity', col = 'black', show.legend = FALSE) +
  theme_bw(base_size = 14) +
  labs(y = 'Number of samples',
       x = 'habitat type') +
  facet_wrap(~method) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#--------------------------------#
# try hierarchical clustering ####
#--------------------------------#

hier_clust <- agnes(ps_wunifrac, method = 'ward')
hier_clust2 <- agnes(ps_wunifrac, method = 'average')
hier_clust3 <- agnes(ps_wunifrac, method = 'complete')
hier_clust4 <- agnes(ps_wunifrac, method = 'single')

hier_clust$ac
hier_clust2$ac
hier_clust3$ac
hier_clust4$ac
# ward method gives the strongest clusters based on agglomerative coefficienct

cor(cophenetic(hier_clust), ps_wunifrac)
cor(cophenetic(hier_clust2), ps_wunifrac)
cor(cophenetic(hier_clust3), ps_wunifrac)
cor(cophenetic(hier_clust4), ps_wunifrac)
# average method gives the best representation of the data

# calculate the optimum number of clusters as earlier
hier_clusgap <- clusGap(d_pcoa_correct, FUN = factoextra::hcut, K.max = 12, B = 100, hc_func = 'agnes', hc_method = 'ward')
hier_clusgap2 <- clusGap(d_pcoa_correct, FUN = factoextra::hcut, K.max = 12, B = 100, hc_method = 'average')
hier_nbclust <- NbClust::NbClust(d_pcoa_correct, method = 'ward.D', max.nc = 12)
summary(hier_nbclust)

d_gap <- data.frame(hier_clusgap$Tab, k=1:nrow(hier_clusgap$Tab)) %>%
  data.frame()
d_gap2 <- data.frame(hier_clusgap2$Tab, k=1:nrow(hier_clusgap2$Tab)) %>%
  data.frame()

# calculate the optimal number of clusters
maxSE(d_gap$gap, d_gap$SE.sim, method = 'firstSEmax')
maxSE(d_gap$gap, d_gap$SE.sim, method = 'Tibs2001SEmax')
maxSE(d_gap2$gap, d_gap2$SE.sim, method = 'firstSEmax')
maxSE(d_gap2$gap, d_gap2$SE.sim, method = 'Tibs2001SEmax')

ggplot(d_gap2, aes(k, gap)) +
  geom_line() +
  geom_point(size = 5) +
  geom_linerange(aes(ymin = gap - SE.sim, ymax = gap + SE.sim)) +
  theme_bw(base_size = 14) +
  labs(y = 'Gap statistic',
       x = 'Number of clusters') +
  scale_x_continuous(breaks = 1:12)

factoextra::fviz_nbclust(d_pcoa_correct, FUN = hcut, method = "wss", hc_method = 'average')
factoextra::fviz_nbclust(d_pcoa_correct, FUN = hcut, method = "wss", hc_method = 'ward.D2')
fviz_nbclust(d_pcoa_correct, FUN = hcut, method = "silhouette", hc_method = 'ward.D2')
fviz_nbclust(d_pcoa_correct, FUN = hcut, method = "silhouette", hc_method = 'average')
fviz_nbclust(d_pcoa_correct, FUN = hcut, method = "gap", hc_method = 'ward.D', maxSE = list(method = 'Tibs2001SEmax', SE.factor = 1))
fviz_nbclust(d_pcoa_correct, FUN = hcut, method = "gap", hc_method = 'average', maxSE = list(method = 'Tibs2001SEmax', SE.factor = 1))

# compare dendograms using dendextend
dend_ward <- as.dendrogram(hier_clust)
dend_average <-as.dendrogram(hier_clust2)
dend_single <- as.dendrogram(hier_clust4)
dend_complete <- as.dendrogram(hier_clust3)
list_of_dends <- dendlist(dend_ward, dend_average, dend_single, dend_complete)
names(list_of_dends) <- c('ward', 'average', 'single', 'complete')
cor.dendlist(list_of_dends)
cor.dendlist(list_of_dends, method = "common")
# ward, average and complete are all pretty similar
entanglement(dend_ward, dend_average)
entanglement(dend_ward, dend_complete)
entanglement(dend_ward, dend_single)

# change labels of all the samples
labels <- rownames_to_column(d_samp, 'sample') %>% 
  select(sample, habitat_group) %>%
  group_by(habitat_group) %>%
  mutate(hab_group = paste(habitat_group, 1:n(), sep = '')) %>%
  ungroup() %>%
  mutate(order = 1:n()) %>%
  left_join(., cols)

labels(dend_ward) <- labels$hab_group[order.dendrogram(dend_ward)]
labels(dend_average) <- labels$hab_group[order.dendrogram(dend_average)]
labels(dend_complete) <- labels$hab_group[order.dendrogram(dend_complete)]
labels(dend_single) <- labels$hab_group[order.dendrogram(dend_single)]
labels_colors(dend_ward) <- labels$col[order.dendrogram(dend_ward)]
labels_colors(dend_average) <- labels$col[order.dendrogram(dend_average)]

# plot a comparison
tanglegram(dend_ward, dend_average, k_labels = 4, k_branches = 4, margin_inner = 12)
tanglegram(dend_ward, dend_average, margin_inner = 12, common_subtrees_color_branches = TRUE,common_subtrees_color_lines = FALSE)


tanglegram(dend_ward, dend_average, k_labels = 3, k_branches = 3, margin_inner = 12)
tanglegram(dend_complete, dend_ward, k_labels = 4, k_branches = 4, margin_inner = 12)

d_samp2 <- rownames_to_column(d_samp, 'label')
plot(dend_average)

# plot the dendogram using ggtree
p1 <- ggtree(hier_clust) %<+% d_samp2 +
  geom_tippoint(aes(color = habitat_group)) +
  geom_tiplab(aes(label = habitat_group), offset = 0.01, size = MicrobioUoE::pts(8)) +
  scale_color_manual(values = setNames(cols$col, cols$group)) +
  ggplot2::xlim(0, 1.5) +
  labs(title = 'Ward.D2 agglomeration')

p2 <- ggtree(hier_clust2) %<+% d_samp2 +
  geom_tippoint(aes(color = habitat_group)) +
  geom_tiplab(aes(label = habitat_group), offset = 0.01, size = MicrobioUoE::pts(8)) +
  scale_color_manual(values = setNames(cols$col, cols$group)) +
  ggplot2::xlim(0, 0.4) +
  labs(title = 'average agglomeration method')

p3 <- ggtree(hier_clust3) %<+% d_samp2 +
  geom_tippoint(aes(color = habitat_group)) +
  geom_tiplab(aes(label = habitat_group), offset = 0.01, size = MicrobioUoE::pts(8)) +
  scale_color_manual(values = setNames(cols$col, cols$group)) +
  ggplot2::xlim(0, 0.5) +
  labs(title = 'complete agglomeration method')

p4 <- ggtree(hier_clust4) %<+% d_samp2 +
  geom_tippoint(aes(color = habitat_group)) +
  geom_tiplab(aes(label = habitat_group), offset = 0.01, size = MicrobioUoE::pts(8)) +
  scale_color_manual(values = setNames(cols$col, cols$group)) +
  ggplot2::xlim(0, 0.25) +
  labs(title = 'single agglomeration method')

p1 + p2 + p3 + p4 + plot_layout(guides = 'collect', ncol = 2) & guides(color = guide_legend(override.aes = list(size = 6)))

# assign clusters for hierarchical clusterings
# see the cluster
clustering <- cutree(hier_clust, k = 9)
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
  #ggforce::geom_mark_hull(aes(PCoA1, PCoA2, col = cluster), d_clustering, show.legend = FALSE, concavity = 2) +
  #ggrepel::geom_label_repel(aes(PCoA1, PCoA2, col = habitat, label = id), data = to_label) +
  geom_point(aes(PCoA1, PCoA2, col = cluster), d_centroids, size = 3, show.legend = FALSE) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = sample, col = cluster), d_lines, show.legend = FALSE) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12)) +
  facet_wrap(~method)

# create dataframe for cluster assignments for each sample
d_clusters <- sc_clustering(d_pcoa_correct, 5)



# try and use scclust 
x <- cmdscale(ps_wunifrac, 1, eig=T)




# Plot the eigenvalues and choose the correct number of dimensions (eigenvalues close to 0)
plot(x$eig, 
     type="h", lwd=5, las=1, 
     xlab="Number of dimensions", 
     ylab="Eigenvalues")

x <- cmdscale(ps_wunifrac, 45)

clustering2 <- pvclust::pvclust(t(x), method.hclust = 'ward.D2')

plot(clustering2)
pvrect(clustering2, alpha = 0.5)
print(clustering2, digits=3)


