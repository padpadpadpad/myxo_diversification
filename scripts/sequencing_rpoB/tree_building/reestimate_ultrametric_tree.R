# make a single tree ultrametric by re-estimating with a single clock model

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
library(phangorn)

# read in tree
# use either just rooted or after adding tiny tips to the branches
tree <- ape::read.nexus('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv.raxml_rooted.tre')
  
ggtree(tree) + geom_tiplab()
  
# check if tree is rooted
is.rooted(tree)
is.ultrametric(tree)

# if minimum edge length is 0, add a tiny number onto this
if(min(tree$edge.length == 0)){
  tree$edge.length[tree$edge.length == min(tree$edge.length)] <- min(tree$edge.length) + 1e-05
}

# look at edge lengths
hist(tree$edge.length)

# remove family names from tip labels
# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

tree$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

# re-estimate ultrametric phylogenetic tree from the topology of the original tree and the sequencing alignment

# read in alignment for myxo 91%
# has to be phyDat object
alignment <- read.FASTA('data/sequencing_rpoB/raxml/alignment/alignment_asv.fasta') %>%
  as.phyDat()

# create distance matrix
# use hamming distance but this could be changes
dist_matrix <- dist.hamming(alignment)

# make tree ultrametric
tree_ultra <- nnls.phylo(tree, dist_matrix, method = 'ultrametric')

ggtree(tree_ultra) + geom_tiplab()

# optimise ultrametric tree using ml
fit <- pml(tree_ultra, alignment, k=4, model = 'GTR')

ggtree(fit$tree) + geom_tiplab()

# optimise ultrametric tree using another method
fit2 <- optim.pml(fit, optRooted = TRUE, model = 'GTR', optGamma = TRUE)

ggtree(fit2$tree) + geom_tiplab()

# save out ultrametric tree
write.tree(fit2$tree, 'data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_phangorn.tre')
