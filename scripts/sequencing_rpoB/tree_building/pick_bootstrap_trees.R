# load in bootstrapped trees, pick 9 random trees and resave them

# load in packages
librarian::shelf(ape, tidyverse)

# set seed
set.seed(42)

# load in trees
trees <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/boots.raxml.bootstraps')

# pick 9 random trees
trees <- trees[sample(1:length(trees), 9)]

# load in treepl tree with family name tip labels
tree1 <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_treepl_cv.tre')

# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

# extract original tip labels
tree1_tiplabels <- tibble(tree1_label = tree1$tip.label,
                          tree2_label = purrr::map_chr(tree1_label, strsplit_mod
                          ))

trees2 <- trees

for(i in 1:9){
  temp <- tibble(tree2_label = trees[[i]]$tip.label) %>%
    left_join(., tree1_tiplabels)
  trees[[i]]$tip.label <- temp$tree1_label
  
  # save out tree
  write.tree(trees[[i]], paste0('data/sequencing_rpoB/raxml/trees/myxo_asv/raxml_bootstrap_', i, '.tree'))
  }

trees[[1]]$tip.label[1:6]
trees2[[1]]$tip.label[1:6]
trees[[9]]$tip.label[1:6]
trees2[[9]]$tip.label[1:6]

