# compare the phylogenetic ASV trees

librarian::shelf(dendextend, ape, TreeDist, tidyverse)

# load in trees
tree1 <- read.tree("data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_treepl_cv_node_labels.tre")
tree2 <- read.tree("data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_1_treepl_cv.tre")
tree3 <- read.tree("data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_2_treepl_cv.tre")
tree4 <- read.tree("data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_3_treepl_cv.tre")
tree5 <- read.tree("data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_4_treepl_cv.tre")
tree6 <- read.tree("data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_5_treepl_cv.tre")
tree7 <- read.tree("data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_6_treepl_cv.tre")
tree8 <- read.tree("data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_7_treepl_cv.tre")
tree9 <- read.tree("data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_8_treepl_cv.tre")
tree10 <- read.tree("data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_9_treepl_cv.tre")

# combine into a dend list
trees <- dendlist(as.dendrogram(tree1), as.dendrogram(tree2), as.dendrogram(tree3), as.dendrogram(tree4), as.dendrogram(tree5), as.dendrogram(tree6), as.dendrogram(tree7), as.dendrogram(tree8), as.dendrogram(tree9), as.dendrogram(tree10))

# look at correlation between trees
d_cophenetic <- cor.dendlist(trees, method = 'cophenetic')
d_cophenetic <- as.data.frame(d_cophenetic) %>%
  mutate(tree = c('best', paste('boot', 1:9, sep = ' ')))
colnames(d_cophenetic)[grepl('V', colnames(d_cophenetic))] <- c('best', paste('boot', 1:9, sep = ' '))

d_cophenetic <- pivot_longer(d_cophenetic, -tree, names_to = 'tree2', values_to = 'correlation') %>%
  filter(tree != tree2) %>%
  group_by(tree) %>%
  summarise(mean = mean(correlation),
            min = min(correlation),
            max = max(correlation))