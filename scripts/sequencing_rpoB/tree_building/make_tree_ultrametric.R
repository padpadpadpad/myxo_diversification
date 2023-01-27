# make trees ultrametric

#--------------------------#
# what this script does ####
#--------------------------#

# for each rerooted tree data/sequencing_rpoB/phyloseq/myxococcus/clustered
# 1. reads in tree
# 2. checks it is rooted
# 3. makes the tree ultrametric
# 4. checks the tree is ultrametric
# 5. saves out the ultrametric tree twice - one with family names in the tip labels and one without

# load in packages
library(phytools)
library(ape)
library(tidyverse)
library(here)

here::i_am('scripts/sequencing_rpoB/tree_building/make_tree_ultrametric.R')

# set percent similarity - those used in asvs_to_otus.R
percent_similarity <- c(99:90, 97.7, 85, 80, 'asv')

# read in tree
# use either just rooted or after adding tiny tips to the branches
tree <- read.tree(here('data/sequencing_rpoB/raxml/trees/myxo_91percent/myxo_91.raxml.reroot2'))
tree <- read.tree(here('data/sequencing_rpoB/raxml/trees/myxo_91percent/myxo_91.raxml.reroot'))

# check if tree is rooted
is.rooted(tree)
is.ultrametric(tree)

# make tree ultrametric using ape::chronos
tree_ultra <- chronos(tree, lambda = 10)
tree_ultra <- chronopl(tree, lambda = 10)



# make tree ultrametric
tree_ultrametric <- force.ultrametric(tree)

is.rooted(tree_ultrametric)


# this does not look great LOL.


#------------------------------#
# old code do not run this! ####
#------------------------------#

# run perl script to fix nonzero branches
# tree_file <- here('data/sequencing_rpoB/raxml/trees/myxo_91percent/myxo_91.raxml.reroot')

# define outfile
# out_file <- here('data/sequencing_rpoB/raxml/trees/myxo_91percent/myxo_91.raxml.reroot2')

# define script file
# script_file <- here('scripts/sequencing_rpoB/tree_building/fixzerobranches.pl')

# run perl script to create constraint tree
# system(paste('perl' , script_file, '-i', tree_file, '>', out_file, sep = ' '))

# tree <- read.tree(here('data/sequencing_rpoB/raxml/trees/myxo_91percent/myxo_91.raxml.reroot2'))
