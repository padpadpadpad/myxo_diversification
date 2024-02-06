# write script to loop through treepl trees and force ultrametric

# load libraries
librarian::shelf(ape, phytools, tidyverse)

# list trees
boot_trees <- list.files('data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/', pattern = 'treepl_cv.tre', full.names = TRUE)

# keep files with the extension .tre at the end of the file name
boot_trees <- boot_trees[grepl('.tre$', boot_trees)]

# loop through trees, force them to be ultrametric, then resave out as the same name
for (i in 1:length(boot_trees)) {
  # read in tree
  tree <- read.tree(boot_trees[i])
  
  # force ultrametric
  tree <- force.ultrametric(tree, message = FALSE)
  
  #cat(is.ultrametric(tree))
  
  # write out tree
  write.tree(tree, boot_trees[i])
}
