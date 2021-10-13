#----------------------------------------#
# Script 2 of 16s sequencing analysis ####
#----------------------------------------#

# what does this script do

# 1. remove amplicons that are >250bp
# 2. looks at read distribution of sequencing run post prevalence filtering
# 3. runs rarefaction curves
# 4. rarefies data
# 5. sets root for tree made in FastTree, chooses the same root for both rarefied and original data

rm(list = ls())

# load packages ####
library(phyloseq)
library(vegan)
library(MicrobioUoE)
library(ggvegan)
library(tidyverse)
library(janitor)
library(doParallel)

# load in extra functions
source('scripts/extra_functions.R')

# set seed
set.seed(42)

# figure path
path_fig <- 'plots/sequencing_16s'

# load data - latest run which we are happy with ####
ps <- readRDS('data/sequencing_16s/ps_16s_prev_filtered.rds')

# show available ranks in the dataset
rank_names(ps)

# fix ranks - remove Species and only keep Species.1 - as per Ben Callahan's recommendation
tax_table(ps) <- tax_table(ps)[,c(1:6, 8)]
colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# get dataframe from phyloseq
d_ps <- speedyseq::psmelt(ps) %>%
  janitor::clean_names() %>%
  group_by(sample) %>%
  mutate(rel_abundance = abundance/sum(abundance)) %>%
  ungroup()

#---------------------------------------#
# 1. look at ASV length distribution ####
#---------------------------------------#

# look at distribution of OTU length
ggplot(d_ps, aes(nchar(otu))) +
  geom_histogram()

summary(nchar(d_ps$otu))

# remova OTUs longer than 250 bps
all_taxa <- taxa_names(ps)
to_drop <- all_taxa[nchar(all_taxa) > 250]
to_keep <- all_taxa[!all_taxa %in% to_drop]

ps <- prune_taxa(to_keep, ps)

d_ps <- speedyseq::psmelt(ps) %>%
  janitor::clean_names() %>%
  group_by(sample) %>%
  mutate(rel_abundance = abundance/sum(abundance)) %>%
  ungroup()

# look at distribution of OTU length
ggplot(d_ps, aes(nchar(otu))) +
  geom_histogram()

#----------------------------------------------------------------#
# 2. look at read distribution following prevalence filtering ####
#----------------------------------------------------------------#

num_reads <- tibble(reads = sample_sums(ps), sample = names(sample_sums(ps)), x = 'blah')

# remove samples with too few reads
label_reads <- filter(num_reads, reads >= quantile(reads, 0.95) | reads <= quantile(reads, 0.05))

pos = position_jitter(width = 0.2, seed = 1)

ggplot(num_reads, aes(x, reads)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(position = position_jitter(width = 0.2), size = 3, fill = 'white', shape = 21, data = filter(num_reads, ! sample %in% label_reads$sample)) +
  geom_point(position = position_jitter(width = 0.2, seed = 1), size = 3, fill = 'white', shape = 21, data = label_reads) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(10^0, 10^6),
                minor_breaks = NULL) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_blank()) +
  labs(x = '',
       y = 'Number of reads',
       title = 'Number of reads per sample') +
  ggrepel::geom_label_repel(aes(label = sample), data = label_reads, position = position_jitter(width = 0.2, seed = 1))

ggsave(file.path(path_fig, 'read_distribution.pdf'), last_plot(), width = 5, height = 9)

#------------------------------#
# 3. run rarefaction curves ####
#------------------------------#

# code taken from https://github.com/joey711/phyloseq/issues/143

# Setting up and registering the cluster
cl <- makeCluster(detectCores(all.tests=TRUE))
registerDoParallel(cl)

to_do <- seq(1, max(num_reads$reads), length.out = 75)

rarefaction_curves <- calculate_rarefaction_curves(ps, c('Observed', 'Shannon'), rep(to_do, each = 10), parallel=TRUE)

stopCluster(cl)

rarefaction_curves <- group_by(rarefaction_curves, Depth, Sample, Measure) %>%
  summarise(diversity_mean = mean(Alpha_diversity),
            sd = sd(Alpha_diversity),
            .groups = 'drop') %>%
  janitor::clean_names() %>%
  mutate(across(where(is.factor), as.character))

meta_data = sample_data(ps) %>% 
  data.frame() %>%
  mutate(sample = paste('sample_s', id, sep = ''))

rarefaction_curves <- left_join(rarefaction_curves, select(meta_data, sample, habitat_group, location))

ggplot(filter(rarefaction_curves, measure == 'Observed')) +
  geom_line(aes(depth, diversity_mean, group = sample)) +
  geom_ribbon(aes(x = depth, ymin = diversity_mean - sd, ymax = diversity_mean + sd, group = sample), alpha = 0.3) +
  facet_wrap(~habitat_group) +
  theme_bw()

ggsave(file.path(path_fig, 'rarefaction_curves_16s.pdf'), last_plot(), width = 10, height = 7)

#-------------------#
# 4. rarefy data ####
#-------------------#

# rarefy these samples down to the lowest number of reads 
ps_rarefy <- rarefy_even_depth(ps)
# save out datasets

#--------------------------------------------#
# 5. assign new root to phylogenetic tree ####
#--------------------------------------------#

# assign reproducible tree root to the dataset
# https://john-quensen.com/r/unifrac-and-tree-roots/
# Sebastian Schmidt proposed choosing an out group based on the longest tree branch terminating in a tip
# https://github.com/joey711/phyloseq/issues/597

# function for picking tree root:
pick_new_outgroup <- function(tree.unrooted){
  # tablify parts of tree that we need.
  treeDT <-
    cbind(
      data.table::data.table(tree.unrooted$edge),
      data.table::data.table(length = tree.unrooted$edge.length)
    )[1:ape::Ntip(tree.unrooted)] %>%ps_1
    cbind(data.table::data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)}

# set root for rarefied sample then the same one for the none-rarefied sample
out_group <- phy_tree(ps_rarefy) %>% pick_new_outgroup()

phy_tree(ps_rarefy) <- ape::root(phy_tree(ps_rarefy), outgroup=out_group, resolve.root=TRUE)
phy_tree(ps) <- ape::root(phy_tree(ps), outgroup=out_group, resolve.root=TRUE)

# save out datasets
saveRDS(ps, 'data/sequencing_16s/ps_16s_low_depth_removed.rds')
saveRDS(ps_rarefy, 'data/sequencing_16s/ps_16s_rarefied.rds')

