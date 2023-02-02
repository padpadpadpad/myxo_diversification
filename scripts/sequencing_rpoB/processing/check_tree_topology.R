# check tree topology for every OTU cut-off

# load packages
library(here)
library(tidyverse)
library(ggtree)
library(ggnewscale)
library(RColorBrewer)
library(patchwork)
library(ape)
library(phytools)
library(MetBrewer)

# set where I am in the project
here::i_am('scripts/sequencing_rpoB/processing/check_tree_topology.R')

# set habitat colours
cols_hab <- met.brewer('Austria', n = 7)
names(cols_hab) <- c('marine_mud', 'freshwater', 'terrestrial', 'freshwater:terrestrial', 'generalist', 'marine_mud:terrestrial', 'freshwater:marine_mud')
hab_labels <- c('marine mud', 'freshwater', 'terrestrial', 'freshwater + terrestrial', 'generalist', 'marine mud + terrestrial', 'freshwater + marine mud')

# constrained families
constrained_families <- c('Myxococcaceae', 'Vulgatibacteraceae', 'Anaeromyxobacteraceae', 'Polyangiaceae', 'Sandaracinaceae', 'Nannocystaceae', 'Haliangiaceae')

# set variable for ASV or otu similarity
otu_similarity <- c(99:91, 97.7, 'asv') %>% sort()

pdf(here('plots/sequencing_rpoB/analyses/tree_topology_check.pdf'), height = 10, width = 5)

for(i in 1:length(otu_similarity)){
  
  # set otu similarity
  temp_otu_similarity = otu_similarity[i]
  temp_otu_similarity2 = paste(otu_similarity[i], 'percent', sep = '')
  if(otu_similarity[i] == 'asv'){temp_otu_similarity2 = 'asv'}
  
  # read in dataset about habitat preferences
  if(temp_otu_similarity == 'asv'){d_habpref <- read.csv(here(paste('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_', temp_otu_similarity, '.csv', sep = '')))}
  
  if(temp_otu_similarity != 'asv'){d_habpref <- read.csv(here(paste('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_', temp_otu_similarity, 'percent.csv', sep = '')))}
  
  # read in rooted tree - tree has been rooted in figtree
  tree <- ape::read.tree(paste('data/sequencing_rpoB/raxml/trees/myxo_', temp_otu_similarity,'/myxo_', temp_otu_similarity, '_chronopl10.tre', sep = ''))
  

  # read in phyloseq object and grab tax table
  d_taxa <- readRDS(here(paste('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_', temp_otu_similarity2,  '_filt.rds', sep = ''))) %>%
    phyloseq::tax_table() %>%
    data.frame() %>%
    janitor::clean_names() %>%
    rownames_to_column('otu')
  
  # create d_meta
  d_meta <- left_join(dplyr::select(d_habpref, otu, habitat_preference, num_present), dplyr::select(d_taxa, otu:family)) %>%
    mutate(family2 = ifelse(is.na(family), 'uncertain', family))
  
  
  # find 9 most common families based on number of tips assigned to them
  d_common <- d_meta %>%
    group_by(family2) %>%
    tally() %>%
    ungroup() %>%
    mutate(., prop = n / sum(n)) %>%
    filter(., family2 %in% constrained_families)
  
  sum(d_common$prop)
  
  # add in column for colouring branches
  d_meta <- mutate(d_meta, family_common = ifelse(family2 %in% c(constrained_families, 'uncertain'), family2, 'unconstrained')) %>%
    rename(tip_label = otu)
  
  # group tip labels together in terms of their order
  to_group <- split(d_meta$tip_label, d_meta$family_common)
  tree2 <- groupOTU(tree, to_group)
  
  # find the mrca of each of the clades
  d_meta2 <- filter(d_meta, family_common %in% constrained_families) %>%
    dplyr::select(family, tip_label) %>%
    group_by(family) %>%
    nest() %>%
    mutate(mrca = NA) %>%
    ungroup()
  
  #for(i in 1:nrow(d_meta2)){
    #d_meta2$mrca[i] <- findMRCA(tree2, tips = d_meta2$data[[i]]$tip_label)
  #}
  
  #d_meta2 <- dplyr::select(d_meta2, family2 = family, mrca) %>%
    #mutate(blank_label = '')
  
  # create colour palettes
  
  # for different families - add in black and grey for uncertain and uncommon respectively
  cols <- c(colorRampPalette(brewer.pal(11, "Spectral"))(nrow(d_common)), 'black', 'grey')
  names(cols) <- c(sort(d_common$family2), 'uncertain', 'unconstrained')
  
  # make a separate column for the rare states to make their size bigger!
  d_meta <- mutate(d_meta, rare = ifelse(habitat_preference %in% c('generalist', 'marine_mud:terrestrial'), 'rare', 'common'))
  
  # first only plot taxonomy, also make non circular so we can see groupings properly
  tree_plot <- ggtree(tree2, aes(col = group)) %<+% filter(d_meta) +
    scale_color_manual('Family (branch colours)', values = cols) +
    guides(color = guide_legend(override.aes = list(linewidth = 3))) +
    #geom_cladelab(data = d_meta2,
                  #mapping = aes(node = mrca,
                                #color = family2,
                                #label = blank_label),
                  #offset = 0.1,
                  #barsize = 2,
                  #show.legend = FALSE) +
    ggtitle(paste('OTU similarity = ', temp_otu_similarity2, sep = ''))
  
  print(tree_plot)
  
}

dev.off()
