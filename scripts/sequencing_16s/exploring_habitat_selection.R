# look at clustering of otus and samples ####

# clean workspace
rm(list = ls())

# load packages ####
library(microViz)
library(phyloseq)
library(vegan)
library(patchwork)
library(tidyverse)
library(palettetown)
library(MicrobioUoE)

# if not installed, install mctoolsr run remotes::install_github('leffj/mctoolsr')

path_fig <- 'sequencing_16S/plots/analyses'

# source extra functions
source('sequencing_16S/scripts/extra_functions.R')

# load data
#ps <- readRDS('data/sequencing/output/run_4/ps_rarefied.rds')
ps <- readRDS('sequencing_16S/data/output/run_merged_runs_new/ps_low_depth_removed.rds')

# replace with new metadata
meta <- read.csv('sequencing_16S/data/metadata_complete_fixed.csv', stringsAsFactors = FALSE)

# read in clustering data


row.names(meta) <- paste('sample_s', meta$id, sep = '')
sample_data(ps) <- sample_data(meta)

summary <- group_by(meta, habitat_group2) %>%
  tally()

sort(sample_sums(ps))
rank_names(ps)

d_ps <- speedyseq::psmelt(ps) %>%
  janitor::clean_names()
d_ps <- group_by(d_ps, sample) %>%
  mutate(total_reads = sum(abundance)) %>%
  ungroup()
d_myxo <- filter(d_ps, phylum == 'Myxococcota') %>%
  group_by(sample) %>%
  mutate(prop = abundance/total_reads) %>%
  ungroup() %>%
  filter(abundance > 0) %>%
  arrange(otu) %>%
  group_by(otu) %>%
  mutate(otu2 = paste('otu', cur_group_id(), sep = '_')) %>%
  ungroup()

otu_1 <- unique(d_myxo$otu[1])

d_myxo_test <- filter(d_myxo, otu == otu_1)

head(d_myxo_test)




d_myxo_summary <- group_by(d_myxo, sample, habitat_group) %>%
  summarise(prop = sum(prop),
            num = n(),
            .groups = 'drop')

d_myxo_summary2 <- group_by(d_myxo, otu) %>%
  summarise(prop = sum(prop),
            num = n(),
            .groups = 'drop')



ggplot(d_myxo_summary, aes(habitat_group, num)) +
  geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(shape = 21, fill = 'white')

ggplot(d_myxo_summary2, aes(prop)) +
  geom_histogram(fill = 'white', col = 'black')
ggplot(d_myxo_summary2, aes(num)) +
  geom_histogram(fill = 'white', col = 'black')

ggplot(d_myxo, aes(habitat, prop)) +
  facet_wrap(~otu2) +
  geom_point()

ps_myxo <- subset_taxa(ps, Phylum == 'Myxococcota')

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
correct_eigenvalues[1:2]

# set up colours based on habitats
cols <- tibble(group = c("woodland_oak", "estuarine mud_low polyhaline", "woodland_pine", "reservoir", "river", "estuarine mud_oligohaline", "estuarine mud_full saline", "beach_supratidal", "pasture", "beach_subtidal","thrift_rhizosphere","beach_seaweed","field_wheat","rock_samphire","marine mud_full saline"),
               col = c('#089a2d', '#995a08', '#106c12', '#1170bd', '#9dcdf4', '#663c05', '#b6966b', '#f5e279', '#61dd1e', '#f2f426', '#c5f8ae', '#a6ab52', '#9ff121', '#5e8128', '#714a03'))
cols <- filter(cols, group != 'thrift_rhizosphere')

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
  scale_color_manual(values = setNames(cols$col, cols$group))

ggsave(file.path(path_fig, 'ordination_16S.pdf'), p1, width = 10, height = 7)
ggsave(file.path(path_fig, 'ordination_16S.png'), p1, width = 10, height = 7)

# look at pairwise differences between groups
# pairwise permanovas
mods <- calc_pairwise_permanovas(ps_wunifrac, d_samp, 'habitat_group_fac', n_perm = 9999)

write.csv(select(mods, X1, X2, R2, pval, pvalFDR), 'sequencing_16S/data/analysis/pairwise_permanovas.csv', row.names = FALSE)
filter(mods, pvalFDR >= 0.05) %>%
  select(., X1, X2, R2, pval, pvalFDR) %>%
  write.csv(., 'sequencing_16S/data/analysis/pairwise_permanovas_notsignificant.csv', row.names = FALSE)

# example for pairwise clusterings
example <- select(mods, X1, X2, pvalFDR) %>%
  mutate(contrast = 1:n(),
         contrast2 = paste(X1, '\nvs.\n', X2, sep = ''),
         contrast2 = gsub('_', ' ', contrast2)) %>%
  pivot_longer(c(X1, X2), names_to = 'drop', values_to = 'group') %>%
  select(-drop)

example_eigenvector <- left_join(example, d_fig$eigenvector) %>%
  filter(contrast <=12)
example_centroids <- left_join(example, d_fig$centroids) %>%
  filter(contrast <=12)
example_lines <- left_join(example, betadisper_lines) %>%
  filter(contrast <=12)

ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = group, shape = location),example_eigenvector) +
  #ggrepel::geom_label_repel(aes(PCoA1, PCoA2, col = habitat, label = id), data = to_label) +
  geom_point(aes(PCoA1, PCoA2, col = group), example_centroids, size = 3) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = sample, col = group), example_lines) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0.5, size = 12),
        legend.position = 'right') +
  scale_color_manual(values = setNames(cols$col, cols$group)) +
  facet_wrap(~contrast2) +
  gghighlight::gghighlight(unhighlighted_params = list(colour = NULL, alpha = 0.1), label_key = habitat_group)

ggsave(file.path(path_fig, 'example_permanovas.pdf'), last_plot(), width = 14, height = 10)
ggsave(file.path(path_fig, 'example_permanovas.png'), last_plot(), width = 14, height = 10)
