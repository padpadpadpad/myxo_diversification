# align sequences for each clustered phyloseq object #
# create a backbone tree for each tree #

#--------------------------#
# what this script does ####
#--------------------------#

# for each phyloseq object in data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered
# 1. reads in phyloseq object
# 2. aligns these sequences
# 3. creates a backbone tree

# load in packages
library(DECIPHER)
library(phyloseq)
library(tidyverse)
library(ggtree)
library(RColorBrewer)
library(ape)

# set percent similarity - those used in asvs_to_otus.R
percent_similarity <- c('asv')

# read in sequence to otu name converter
asv_names <- read_csv('data/sequencing_rpoB/phyloseq/myxococcus/myxo_seq_name_conversion.csv')
asv_names

otu_similarity = 'asv'

# set backbone tree
backbone <- 'data/sequencing_rpoB/raxml/backbone.tre'

# load in myxo phyoseq object
ps_myxo <- readRDS(paste('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_', otu_similarity,  '_filt.rds', sep = ''))

# keep only ASVs that are present in this dataset 
df_seqs <- filter(asv_names, otu_name %in% taxa_names(ps_myxo))

# DNA string set
seqs <- DNAStringSet(df_seqs$seq)
names(seqs) <- df_seqs$otu_name
seqs <- OrientNucleotides(seqs)

# build guide tree
guide_tree <- lapply(order(width(seqs), decreasing=TRUE),
                     function(x) {
                       attr(x, "height") <- 0
                       attr(x, "label") <- names(seqs)[x]
                       attr(x, "members") <- 1L
                       attr(x, "leaf") <- TRUE
                       x
                     })

attr(guide_tree, "height") <- 0.5
attr(guide_tree, "members") <- length(seqs)
class(guide_tree) <- "dendrogram"

# align sequences
alignment <- AlignSeqs(seqs, guideTree = guide_tree, anchor = NA, processors = 3)

# try and run the perl script to create a backbone tree using perl

# get taxonomy table from ps object
tax_table <- tax_table(ps_myxo) %>%
  data.frame() %>%
  rownames_to_column(var = 'tip_label') %>%
  janitor::clean_names()

head(tax_table)

# keep only 15 of each of the taxa within each family
tax_table_sub <- mutate(tax_table, family = ifelse(family %in% c('Haliangiaceae', 'Myxococcaceae', 'Vulgatibacteraceae', 'Anaeromyxobacteraceae', 'Polyangiaceae', 'Sandaracinaceae', 'Nannocystaceae', 'Haliangiaceae'), family, 'other')) %>%
  group_by(family) %>%
  sample_n(15) %>%
  ungroup()

# make this a temp file
temp_file <- tempfile()
write.csv(tax_table_sub, temp_file, row.names = FALSE)

# define output file
out_file <- paste('data/sequencing_rpoB/raxml/constraint_trees/constraint_tree_', otu_similarity,  '_subset.tre', sep = '')

# run perl script to create constraint tree
system(paste('perl' ,'~/Desktop/myxo_diversification/scripts/sequencing_rpoB/tree_building/csv2constraint.pl', '-i', temp_file, '-b', backbone, '>', out_file, sep = ' '))

# filter alignment to just contain sequences in the subsetted taxon table
alignment <- alignment[names(alignment) %in% tax_table_sub$tip_label]

# write alignment to file
writeXStringSet(alignment, 'data/sequencing_rpoB/raxml/alignment/alignment_asv_subset.fasta')

# download temp file to convert tree to format for fasttree for making a constraint
system('perl ~/Downloads/TreeToConstraints.pl < data/sequencing_rpoB/raxml/constraint_trees/constraint_tree_asv_subset.tre > data/sequencing_rpoB/raxml/constraint_trees/constraint_tree_asv_subset.fasttree')

# run fasttree to create a very simple tree consistent with constraints
system('fasttree -nt -constraints data/sequencing_rpoB/raxml/constraint_trees/constraint_tree_asv_subset.fasttree < data/sequencing_rpoB/raxml/alignment/alignment_asv_subset.fasta > data/sequencing_rpoB/raxml/constraint_trees/start_tree_asv.tre')

tree <- ape::read.tree('data/sequencing_rpoB/raxml/constraint_trees/start_tree_asv.tre')

# drop any OTUs not in constrained families
drop.tip(tree, tree$tip.label[!tree$tip.label %in% filter(tax_table_sub, family %in% constrained_families)$tip_label]) %>%
  write.tree('data/sequencing_rpoB/raxml/constraint_trees/start_tree_asv_2.tre')

# two OTUs on different sides of the tree
# otu_3541 - Myxococcaceae
# otu_160 - Nannocystaceae

# check where the node is where the split is
ape::getMRCA(tree, c('otu_3541', 'otu_160'))

# re root to that point
tree2 <- ape::root(tree, node = 122)
ape::is.rooted(tree2)

# plot tree with taxonomic info
constrained_families <- c('Myxococcaceae', 'Vulgatibacteraceae', 'Anaeromyxobacteraceae', 'Polyangiaceae', 'Sandaracinaceae', 'Nannocystaceae', 'Haliangiaceae')

d_meta <- tibble(tip_label = tree2$tip.label) %>%
  left_join(tax_table, by = 'tip_label') %>%
  mutate(family = ifelse(family %in% constrained_families, family, 'other'))

# find 9 most common families based on number of tips assigned to them
d_common <- d_meta %>%
  group_by(family) %>%
  tally() %>%
  ungroup() %>%
  mutate(., prop = n / sum(n)) %>%
  filter(., family %in% constrained_families)

# group tip labels together in terms of their order
to_group <- split(d_meta$tip_label, d_meta$family)
tree3 <- groupOTU(tree2, to_group)

# for different families - add in black and grey for uncertain and uncommon respectively
cols <- c(colorRampPalette(brewer.pal(11, "Spectral"))(nrow(d_common)), 'grey')
names(cols) <- c(sort(d_common$family), 'other')


# first only plot taxonomy, also make non circular so we can see groupings properly
tree_plot <- ggtree(tree3, aes(col = group)) %<+% filter(d_meta) +
  scale_color_manual('Family (branch colours)', values = cols) +
  guides(color = guide_legend(override.aes = list(linewidth = 3))) +
  #geom_cladelab(data = d_meta2,
  #mapping = aes(node = mrca,
  #color = family2,
  #label = blank_label),
  #offset = 0.1,
  #barsize = 2,
  #show.legend = FALSE) +
  NULL

tree_plot

ape::write.tree(tree2, 'data/sequencing_rpoB/raxml/constraint_trees/start_tree_asv_rooted.tre')

