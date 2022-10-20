# make tree ultrametric

library(phytools)


# read in tree
tree <- read.tree('data/sequencing_rpoB/treepl/myxo_treepl_smooth10.tree')

# check if ultrametric
is.ultrametric(tree)

# make tree ultrametric
tree_ultrametric <- force.ultrametric(tree)

write.tree(tree_ultrametric, 'data/sequencing_rpoB/bamm/myxo_treepl_smooth10_ultrametric.tree')
