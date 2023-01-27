# add family to otu names to allow grouping in FigTree to root the tree properly

#--------------------------#
# what this script does ####
#--------------------------#

# load in packages
library(phyloseq)
library(tidyverse)
library(here)

here::i_am('scripts/sequencing_rpoB/tree_building/append_tiplabels.R')

# set percent similarity - those used in asvs_to_otus.R
percent_similarity <- c(99:90, 97.7, 85, 80, 'asv')

# read in sequence to otu name converter
asv_names <- read_csv(here('data/sequencing_rpoB/phyloseq/myxococcus/myxo_seq_name_conversion.csv'))
asv_names

# run for loop to align sequences
for(i in 1:length(percent_similarity)){
  
  # define which level of OTU similarity we are using
  otu_similarity <- paste(percent_similarity[i], 'percent', sep = '')
  if(percent_similarity[i] == 'asv'){otu_similarity = 'asv'}
  
  # load in myxo phyoseq object
  ps_myxo <- readRDS(here(paste('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_', otu_similarity,  '_filt.rds', sep = '')))
  
  # get taxonomy table from ps object
  tax_table <- tax_table(ps_myxo) %>%
    data.frame() %>%
    rownames_to_column(var = 'tip_label') %>%
    janitor::clean_names()
  
  # make this a temp file
  temp_file <- tempfile()
  write.csv(tax_table, temp_file, row.names = FALSE)
  
  # define tree file
  in_tree <- here(paste('data/sequencing_rpoB/raxml/trees/myxo_', otu_similarity, '/myxo_', percent_similarity[i],'.raxml.bestTree', sep = ''))
  
  # define output file
  out_tree <- here(paste('data/sequencing_rpoB/raxml/trees/myxo_', otu_similarity, '/myxo_', percent_similarity[i],'.raxml.family', sep = ''))
  
  # define script file
  script_file <- here('scripts/sequencing_rpoB/tree_building/extend_labels.pl')
  
  # run perl script to create constraint tree
  system(paste('perl' , script_file, '-o', temp_file, '-t', in_tree, '>', out_tree, sep = ' '))
  
}




