#----------------------------------------#
# Script 1 of 16s sequencing analysis ####
#----------------------------------------#

# what does this script do

# 1. fixes metadata
# 2. removes any ASVs not assigned to at least the phylum level
# 3. removes low prevalence ASVs - defined as present in <5% samples and a total read abundance <200
# 4. saves out objects and plots some phylogenetic trees

# load packages
library(patchwork)
library(ggplot2)
library(phyloseq)
library(dplyr)
library(tidyr)

# set seed
set.seed(42)

# figure path
path_fig <- 'plots/sequencing_16s'

# load data - latest run which we are happy with ####
ps <- readRDS('data/sequencing_16s/ps_16s_complete.rds')

sample_names(ps)

#----------------------------#
# 1. fix metadata for 16s ####
#----------------------------#

# replace with new metadata
meta_new <- read.csv('data/metadata.csv', stringsAsFactors = FALSE)
row.names(meta_new) <- paste('sample_s', meta_new$id, sep = '')
sample_data(ps) <- sample_data(meta_new)

# check which samples are not present
filter(meta_new, !row.names(meta_new) %in% sample_names(ps))
# samples missing
# 45 - beach_seaweed

# look at object
ps
# SO MANY TREE NODES

sort(sample_sums(ps))

# first...
rank_names(ps)

# look number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)
# lots of NAs here. These are probably artifacts and can be removed

# raw tree
psraw_tree <- plot_tree(ps, method = "treeonly",
                        ladderize = "left",
                        title = "Raw tree")

#-----------------------------------#
# 2. remove poorly assigned ASVs ####
#-----------------------------------#

# remove NA characterisation
ps0 <- subset_taxa(ps, !is.na(Phylum))
table(tax_table(ps0)[, "Phylum"], exclude = NULL) %>% sort()

# plot this tree
ps0_tree <- plot_tree(ps0, method = "treeonly",
          ladderize = "left",
          title = "Tree with NA phyla removed")

# explore prevalence of the dataframe
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(prevalence = prevdf,
                    total_abundance = taxa_sums(ps0),
                    tax_table(ps0), stringsAsFactors = FALSE) %>%
  mutate(., feature = row.names(.))

# head(prevdf)

# Are there phyla that are comprised of mostly low-prevalence features? Compute the total and average prevalences of the features in each phylum.
prev_sum <- group_by(prevdf, Phylum) %>%
  summarise(., mean_prevalence = mean(prevalence)/nsamples(ps0)*100,
            total_prevalence = sum(prevalence)) %>%
  ungroup()

prev_sum

# get unique taxa for prevalence and plot
prevdf1 <- filter(prevdf, Phylum %in% get_taxa_unique(ps0, "Phylum"))

# plot prevalence of samples by 
ggplot(prevdf1, aes(total_abundance, prevalence / nsamples(ps0),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.1, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~ Phylum) + theme(legend.position="none")

# save plot
ggsave(file.path(path_fig, 'prevalence_phyla_plot.pdf'), last_plot(), height = 10, width = 13)

#----------------------------#
# 3. prevalence filtering ####
#----------------------------#

# delete things under 5% prevalence
# min total read abundance = 200
ps_prevfilt <- microViz::tax_filter(ps0, min_prevalence = 0.05,
                          min_total_abundance = 200)

# plot tree again
ps_prevfilt_tree <- plot_tree(ps_prevfilt, method = "treeonly",
                              ladderize = "left",
                              title = "Tree after prevalence filtering")

# agglomerate taxa by Genus ####

# How many genera would be present after filtering?
length(get_taxa_unique(ps_prevfilt, taxonomic.rank = "Genus"))

# agglomerate taxa
ps_genus = tax_glom(ps_prevfilt, "Genus", NArm = TRUE)

# plot tree again
ps_genus_tree <- plot_tree(ps_genus, ladderize="left",
                           base.spacing = 0, 
                           text.size = 2, 
                           col = 'Phylum',
                           title = "Tree after pooling samples by genus")

ps_genus_tree

#----------------------#
# 4. plot all trees ####
#----------------------#

# plot side by side
plot_all <- psraw_tree + ps0_tree + ps_prevfilt_tree + ps_genus_tree + theme(legend.position = 'none') + plot_layout(ncol = 2)
ggsave(file.path(path_fig, 'tree_compare_all.pdf'), plot_all, height = 14, width = 14)

# plot with taxa colours
ggsave(file.path(path_fig, 'tree_phylum.pdf'), ps_genus_tree, height = 14, width = 12)

# save phyloseq object that has been prevalence filtered
saveRDS(ps_prevfilt, 'data/sequencing_16s/ps_16s_prev_filtered.rds')
