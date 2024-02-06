# try and set bamm priors

librarian::shelf(BAMMtools, ape)

# load in tree
tree <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_treepl_cv_node_labels.tre')

setBAMMpriors(tree, outfile = 'data/sequencing_rpoB/bamm/revision/asv/bamm_priors.txt')
