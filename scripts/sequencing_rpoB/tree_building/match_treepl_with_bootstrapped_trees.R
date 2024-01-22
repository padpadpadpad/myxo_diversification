#---------------------------------------------------------------------#
# match nodes between best rooted raxml tree and bootstrapped tree ####
#---------------------------------------------------------------------#

# load in packages
librarian::shelf(phytools, ape, tidyverse)

# load in best raxml tree from treePL CV
tree1 <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_treepl_cv.tre')

# load in bootstrapped tree with tbe values from raxml - unrooted
tree2 <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv.raxml_tbe.tree')

# append family name onto tip label for tree2 ####

# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

# extract original tip labels
tree1_tiplabels <- tibble(tree1_label = tree1$tip.label,
                          tree2_label = purrr::map_chr(tree1_label, strsplit_mod
))

tree_tiplabels <- tibble(tree2_label = tree2$tip.label) %>%
  left_join(., tree1_tiplabels)

tree2$tip.label <- tree_tiplabels$tree1_label

# save out tree
write.tree(tree2, 'data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv.raxml_tbe_v2.tree')

# root tree and save out as NEXUS file in FigTree
# order tips by descending order
# use save as currently displayed and include annotations

# load in tree2 now after rooting correctly
tree2 <- readNexus('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv.raxml_tbe_rooted.tree', format = 'raxml')

# check node labels exist
tree2$node.label

# extract tree data from tree2
# keep node labels only
tr2 <- treeio::as_tibble(tree2) %>%
  filter(!str_detect(label, 'otu')) %>%
  dplyr::select(tree2_node = node, label) %>%
  mutate(tbe = parse_number(label))

# look at distribution of tbe values
hist(tr2$tbe)

# do all the tip names match up between the treePL tree and the bootstrapped tree?
sum(tree1$tip.label == tree2$tip.label) == length(tree1$tip.label)

tip_labels <- tibble(tree1_label = tree1$tip.label,
                     tree2_label = tree2$tip.label)

# no
filter(tip_labels, tree1_label != tree2_label)

# remove the text to make sure all the tip labels match up
tree2$tip.label <- gsub('\\[&!rotate=false]', '', tree2$tip.label)

# recheck it
sum(tree1$tip.label == tree2$tip.label) == length(tree1$tip.label)
# TRUE = success

tree1
tree2

# try and match nodes between the two trees
node_match <- matchNodes(tree1, tree2) %>%
  data.frame()

head(node_match)

# do they ever not match?
filter(node_match, tr1 != tr2)
# No! perfect match

# can just add node labels into tree1
tree1$node.label <- parse_number(tree2$node.label) %>% as.character()

# save out tree and keep node labels
write.tree(tree1, 'data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_treepl_cv_node_labels.tre')
