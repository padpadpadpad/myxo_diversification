# try and set bamm priors

librarian::shelf(BAMMtools, ape)

# load in tree
tree <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_95/myxo_95_treepl_cv_node_labels.tre')

ape::ltt.plot(tree, log = 'y')

setBAMMpriors(tree, outfile = 'data/sequencing_rpoB/bamm/revision/myxo95/bamm_priors.txt')

# load in tree
tree <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_97.7/myxo_97.7_treepl_cv_node_labels.tre')

ape::ltt.plot(tree, log = 'y')

setBAMMpriors(tree, outfile = 'data/sequencing_rpoB/bamm/revision/myxo97.7/bamm_priors.txt')
