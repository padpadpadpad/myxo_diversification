# make trees ultrametric

#--------------------------#
# what this script does ####
#--------------------------#

# for each rerooted tree data/sequencing_rpoB/phyloseq/myxococcus/clustered
# 1. reads in tree
# 2. checks it is rooted
# 3. makes the tree ultrametric by re-estimating the tree using the topology of the original raxml tree and the sequence alignment
# 4. checks the tree is ultrametric
# 5. saves out the ultrametric tree twice - one with family names in the tip labels and one without

# load in packages
library(phytools)
library(phangorn)
library(tidyverse)
library(here)

here::i_am('scripts/sequencing_rpoB/tree_building/make_tree_ultrametric.R')

# set percent similarity - those used in asvs_to_otus.R
percent_similarity <- c(99:90, 97.7, 85, 80, 'asv')

# read in tree
# use either just rooted or after adding tiny tips to the branches
tree <- read.tree(here('data/sequencing_rpoB/raxml/trees/myxo_91percent/myxo_91.raxml.reroot'))
#tree_root_old <- read.tree(here('data/sequencing_rpoB/raxml/rerooted-pruned.tre'))
#tree_ultra_old <- read.tree(here('data/sequencing_rpoB/raxml/rerooted-pruned-chronopl10.tre'))

plot(tree)

# check if tree is rooted
is.rooted(tree)
is.ultrametric(tree)

# re-estimate ultrametric phylogenetic tree from the topology of the original tree and the sequencing alignment

# read in alignment for myxo 91%
# has to be phyDat object
alignment <- read.FASTA(here('data/sequencing_rpoB/raxml/alignment/alignment_91percent.fasta')) %>%
  as.phyDat()

# create distance matrix
# use hamming distance but this could be changes
dist_matrix <- dist.hamming(alignment)

# alter tip labels to remove family as they will not link to the distance matrix
# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))
tree2 <- tree
tree2$tip.label <- purrr::map_chr(tree2$tip.label, strsplit_mod)

# make tree ultrametric
tree_ultra <- nnls.phylo(tree2, dist_matrix, method = 'ultrametric')

# optimise ultrametric tree using ml
fit <- pml(tree_ultra, alignment, k=4, model = 'GTR')
# optimise ultrametric tree using another method
fit2 <- optim.pml(fit, optRooted = TRUE, model = 'GTR', optGamma = TRUE)

# check tree is still rooted
is.rooted(fit2$tree)
# check tree is now ultrametric
is.ultrametric(fit2$tree)

# save out ultrametric tree
write.tree(fit2$tree, here('data/sequencing_rpoB/raxml/trees/myxo_91percent/myxo_91.phangorn.ultrametric.tre'))

#------------------------------#
# old code do not run this! ####
#------------------------------#

# does work for the old tree - albeit it takes a long time
tree_ultra_old <- chronopl(tree_root_old, lambda = 10)

# try make tree smaller to see if it works
tree_sub <- keep.tip(tree, sample(tree$tip.label, 750))

# look at distribution of branch lengths to see if there is a difference between the old tree and the new tree
hist(tree$edge.length)
hist(tree2$edge.length)
hist(tree_root_old$edge.length)
hist(tree_ultra_old$edge.length)
hist(fit2$tree$edge.length)


# now only on the tiny lengths
hist(tree$edge.length[tree$edge.length<0.04])
hist(tree2$edge.length[tree2$edge.length<0.04])
hist(tree_root_old$edge.length[tree_root_old$edge.length<0.04])
hist(tree_ultra_old$edge.length[tree_ultra_old$edge.length<0.04])

# look at the minimum branch lengths
min(tree$edge.length)
min(tree2$edge.length)
min(tree_root_old$edge.length)
min(tree_ultra_old$edge.length)

# make tree ultrametric using ape::chronos
tree_ultra <- chronos(keep.tip(tree, sample(tree$tip.label, 750)), lambda = 10)

plot(tree_ultra)

# tips not in this tree
tips_binned <- tree$tip.label[!tree$tip.label %in% tree_ultra$tip.label]

# plot whether tips was included
d <- data.frame(tip_label = tree$tip.label) %>%
  mutate(included = ifelse(tip_label %in% tips_binned, 'no', 'yes'))

ggtree(tree) %<+% d + 
  geom_tiplab(aes(color=included))

tree_ultra <- chronopl(keep.tip(tree, sample(tree$tip.label, 750)), lambda = 10)

# make tree ultrametric using ape::chronos
# doesnt work for the 91% tree
# tree_ultra <- chronos(tree2, lambda = 10, model = 'correlated')
# tree_ultra <- chronopl(tree, lambda = 10)
# Cannot understand why it does work sometimes but not at others. It is not because of zero length terminal branches as I have checked their distribution and there are not any.

# make tree ultrametric
#tree_ultrametric <- force.ultrametric(tree)

#is.rooted(tree_ultrametric)

# run perl script to fix nonzero branches
# tree_file <- here('data/sequencing_rpoB/raxml/trees/myxo_91percent/myxo_91.raxml.reroot')

# define outfile
# out_file <- here('data/sequencing_rpoB/raxml/trees/myxo_91percent/myxo_91.raxml.reroot2')

# define script file
# script_file <- here('scripts/sequencing_rpoB/tree_building/fixzerobranches.pl')

# run perl script to create constraint tree
# system(paste('perl' , script_file, '-i', tree_file, '>', out_file, sep = ' '))

# tree <- read.tree(here('data/sequencing_rpoB/raxml/trees/myxo_91percent/myxo_91.raxml.reroot2'))
