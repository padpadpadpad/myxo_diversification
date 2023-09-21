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
library(here)

here::i_am('scripts/sequencing_rpoB/processing/align_myxo.R')

# set percent similarity - those used in asvs_to_otus.R
percent_similarity <- c(99:90, 97.7, 85, 80, 'asv')

# read in sequence to otu name converter
asv_names <- read_csv(here('data/sequencing_rpoB/phyloseq/myxococcus/myxo_seq_name_conversion.csv'))
asv_names

# set backbone tree
backbone <- here('data/sequencing_rpoB/raxml/backbone.tre')

# run for loop to align sequences
for(i in 1:length(percent_similarity)){
  
  # define which level of OTU similarity we are using
  otu_similarity <- paste(percent_similarity[i], 'percent', sep = '')
  if(percent_similarity[i] == 'asv'){otu_similarity = 'asv'}
  
  # load in myxo phyoseq object
  ps_myxo <- readRDS(here(paste('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_', otu_similarity,  '_filt.rds', sep = '')))
  
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
  
  # save out alignment
  writeXStringSet(alignment, file=here(paste('data/sequencing_rpoB/raxml/alignment/alignment_', otu_similarity,  '.fasta', sep = '')), compress = FALSE)
  
  # try and run the perl script to create a backbone tree using perl
  
  # get taxonomy table from ps object
  tax_table <- tax_table(ps_myxo) %>%
    data.frame() %>%
    rownames_to_column(var = 'tip_label') %>%
    janitor::clean_names()
  
  # make this a temp file
  temp_file <- tempfile()
  write.csv(tax_table, temp_file, row.names = FALSE)
  
  # define output file
  out_file <- here(paste('data/sequencing_rpoB/raxml/constraint_trees/constraint_tree_', otu_similarity,  '.tre', sep = ''))
  
  # run perl script to create constraint tree
  system(paste('perl' ,'~/Desktop/myxo_diversification/scripts/sequencing_rpoB/tree_building/csv2constraint.pl', '-i', temp_file, '-b', backbone, '>', out_file, sep = ' '))
  
}




