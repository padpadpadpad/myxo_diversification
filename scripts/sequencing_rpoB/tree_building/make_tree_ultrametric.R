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



# re-estimate ultrametric phylogenetic tree from the topology of the original tree and the sequencing alignment

# read in alignment for myxo 91%
# has to be phyDat object
alignment <- read.FASTA(here('data/sequencing_rpoB/raxml/alignment/alignment_91percent.fasta')) %>%
  as.phyDat()

# create distance matrix
# use hamming distance but this could be changes
dist_matrix <- dist.hamming(alignment)


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
