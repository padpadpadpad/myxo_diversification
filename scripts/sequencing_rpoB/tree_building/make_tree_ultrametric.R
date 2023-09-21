# make trees ultrametric

#--------------------------#
# what this script does ####
#--------------------------#

# for each rerooted tree data/sequencing_rpoB/phyloseq/myxococcus/clustered
# 1. reads in tree
# 2. checks it is rooted
# 3. makes the tree ultrametric by using chronopl()
# 4. checks the tree is ultrametric
# 5. saves out the ultrametric tree after removing family names

# load in packages
library(phytools)
library(tidyverse)
library(ggtree)
library(here)

here::i_am('scripts/sequencing_rpoB/tree_building/make_tree_ultrametric.R')

# set percent similarity - those used in asvs_to_otus.R
percent_similarity <- c(99:91, 97.7, 'asv')

# run for loop to align sequences
for(i in 1:length(percent_similarity)){
  
  # read in tree
  # use either just rooted or after adding tiny tips to the branches
  tree <- read.nexus(here(paste('data/sequencing_rpoB/raxml/trees/myxo_', percent_similarity[i], '/myxo_', percent_similarity[i], '.raxml_rooted.tre', sep = '')))
  
  ggtree(tree) + geom_tiplab()
  
  # check if tree is rooted
  is.rooted(tree)
  is.ultrametric(tree)
  
  # if minimum edge length is 0, add a tiny number onto this
  if(min(tree$edge.length == 0)){
    tree$edge.length[tree$edge.length == min(tree$edge.length)] <- min(tree$edge.length) + 1e-05
  }
 
  # try using chronos
  #tree_ultra <- chronos(tree, lambda = 10)
  tree_ultra <- chronopl(tree, lambda = 10)
  
  ggtree(tree_ultra) + geom_tiplab()
  
  # look at edge lengths
  hist(tree$edge.length)
  hist(tree_ultra$edge.length)
  
  is.ultrametric(tree_ultra)
  
  # alter tip labels to remove family as they will not link to the distance matrix
  # write function to remove family labels
  strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))
  
  tree_ultra$tip.label <- purrr::map_chr(tree_ultra$tip.label, strsplit_mod)
  
  # save out ultrametric tree
  write.tree(tree_ultra, here(paste('data/sequencing_rpoB/raxml/trees/myxo_', percent_similarity[i], '/myxo_', percent_similarity[i], '_chronopl10.tre', sep = '')))
  
}