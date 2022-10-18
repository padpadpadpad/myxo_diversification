#-------------------------------------------#
# look at clustering of otus and samples ####
#-------------------------------------------#

# what does this script do

# 1. run PCoA on all samples based on habitat
# 2. run pairwise permanovas on all habitats to see which are not significantly different

# clean workspace
rm(list = ls())

# load packages ####
library(microViz)
library(phyloseq)
library(vegan)
library(patchwork)
library(tidyverse)
library(palettetown)
library(MicrobioUoE) # remotes::install_github('padpadpadpad/MicrobioUoE')
# if not installed, install mctoolsr run remotes::install_github('leffj/mctoolsr')

path_fig <- 'plots/sequencing_16s'

# source extra functions
source('scripts/extra_functions.R')

# load data
ps <- readRDS('data/sequencing_16s/ps_16s_low_depth_removed.rds')

# transform counts to relative abundances for ordination
ps_prop <- transform_sample_counts(ps, function(x){x / sum(x)})

#-----------------------------------------------------------------------#
# 1. Run PCoA looking at community composition of different habitats ####
#-----------------------------------------------------------------------#

# wrangle the metadata
d_samp <- data.frame(sample_data(ps_prop))
d_samp <- mutate(d_samp, location_fac = as.factor(location),
                 habitat_group_fac = as.factor(habitat_group)) %>%
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
betadisper_lines <- left_join(betadisper_lines, rownames_to_column(d_samp, var = 'sample') %>% select(sample, location))
d_fig$eigenvector <- left_join(d_fig$eigenvector, rownames_to_column(d_samp, var = 'sample') %>% select(sample, location)) 

d_fig$eigenvalue <- mutate(d_fig$eigenvalue, percent = eig/sum(eig)*100)

# get correct eigenvalue
correct_eigenvalues <- ape::pcoa(ps_wunifrac, correction = 'cailliez', d_samp$id) %>%
  .$values %>%
  pull(Rel_corr_eig)
correct_eigenvalues[1:2]

# set up colours based on habitats
cols <- tibble(group = c("woodland_oak", "estuarine mud_low polyhaline", "woodland_pine", "reservoir", "river", "estuarine mud_oligohaline", "estuarine mud_full saline", "beach_supratidal", "pasture", "beach_subtidal","thrift_rhizosphere","beach_seaweed","field_wheat","rock_samphire","marine mud_full saline"),
               col = c('#089a2d', '#995a08', '#106c12', '#1170bd', '#9dcdf4', '#663c05', '#b6966b', '#f5e279', '#61dd1e', '#f2f426', '#c5f8ae', '#a6ab52', '#9ff121', '#5e8128', '#714a03'),
               hab_order = c(1.1, 2.1, 1.2, 3.1, 3.2, 2.2, 2.3, 2.4, 1.3, 2.5, 2.6, 2.7, 1.4, 1.5, 2.8))
cols <- cols %>% arrange(hab_order)

# plot PCoA
p1 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = group, shape = location),d_fig$eigenvector) +
  geom_point(aes(PCoA1, PCoA2, col = group), d_fig$centroids, size = 5, show.legend = FALSE) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = sample, col = group), betadisper_lines) +
  #theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12),
        legend.position = 'right') +
  scale_color_manual('Sampled habitat', values = setNames(cols$col, cols$group)) +
  scale_shape_discrete('Sampled location') +
  labs(y = 'Axis 2 (12.78%)',
       x = 'Axis 1 (35.26%)')

ggsave(file.path(path_fig, 'ordination_16S.pdf'), p1, width = 10, height = 7)
ggsave(file.path(path_fig, 'ordination_16S.png'), p1, width = 10, height = 7)

p1 + MicrobioUoE::theme_black(base_size = 24) +
  theme(plot.background = element_rect(fill = '#404040', color = '#404040'),
        panel.background = element_rect(fill = '#404040'),
        legend.position = "none",
        panel.grid.major = element_blank())
  
ggsave(file.path(path_fig, 'ordination_16S_presentation.pdf'), last_plot(), width = 10, height = 7)


#-----------------------------------------------------#
# 2. look at pairwise differences between habitats ####
#-----------------------------------------------------#

mods <- calc_pairwise_permanovas(ps_wunifrac, d_samp, 'habitat_group_fac', n_perm = 9999)

write.csv(select(mods, X1, X2, R2, pval, pvalFDR), 'data/sequencing_16s/pairwise_permanovas.csv', row.names = FALSE)
filter(mods, pvalFDR >= 0.05) %>%
  select(., X1, X2, R2, pval, pvalFDR) %>%
  write.csv(., 'data/sequencing_16s/pairwise_permanovas_notsignificant.csv', row.names = FALSE)

filter(mods, pvalFDR >= 0.05) %>%
  select(., X1, X2, R2, pval, pvalFDR)

# Habitats that were not significantly different
# thrift rhizosphere from almost everything
# beach supratidal from almost everything
# estuarine mud_low polyhaline - estuarine mud oligohaline
# field wheat - pasture