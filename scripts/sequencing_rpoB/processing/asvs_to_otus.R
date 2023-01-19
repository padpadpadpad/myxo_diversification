# Cluster existing phyloseq object to a specific OTU similarity.
# Also do prevalence filtering

#--------------------------#
# what this script does ####
#--------------------------#

# reads in myxo asv object
# aligns this object and creates a disatnace matrix
# clusters the ASVs at a given similarity
# saves out clustered phyloseq object
# does prevalence filtering: present in at least four samples and have a total abundance > 100

# load in packages
library(phyloseq)
library(dada2)
library(DECIPHER)
library(here)
library(tidyverse)

# set number of processors to use
num_processors <- 4

# set where we are
here::i_am('scripts/sequencing_rpoB/processing/asvs_to_otus.R')

# read in phyloseq object
ps <- readRDS(here('data/sequencing_rpoB/phyloseq/myxococcus/clustered/ps_otu_asv.rds'))

# prevalence filter the raw ASV data
ps_sub <- microViz::tax_filter(ps, min_prevalence = 4, min_total_abundance = 100)

# save this out
saveRDS(ps_sub, here('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_asv_filt.rds'))

# read in myxo seq name conversion
seqs_df <- read.csv(here('data/sequencing_rpoB/phyloseq/myxococcus/myxo_seq_name_conversion.csv')) %>%
  filter(otu_name %in% taxa_names(ps))

seqs <- DNAStringSet(seqs_df$seq)
names(seqs) <- seqs_df$otu_name # This propagates to the tip labels of the tree 

# DNA string set
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

# align sequences - this takes a long time on a single machine
alignment <- AlignSeqs(seqs, guideTree = guide_tree, anchor = NA, processors = num_processors)

# calculate distance matrix for each sequence
dist_matrix <- DECIPHER::DistanceMatrix(alignment, processors = num_processors)

# set percent similarity
percent_similarity <- c(99:90, 97.7, 85, 80)
cut_off = (100-percent_similarity) / 100

# run a for loop to cluster the samples for each percent similarity
for(i in 1:length(cut_off))
{
  
  # run cluster algorithm - use the UPGMA algorithm
  # replacement for IdClusters is shown here: https://github.com/benjjneb/dada2/issues/947#issuecomment-1277776614
  clusters <- DECIPHER::TreeLine(
    myDistMatrix = dist_matrix, 
    method = "UPGMA",
    cutoff = cut_off[i],
    type = "clusters",
    processors = num_processors
  )
  
  # create new phyloseq object from the 
  ps0 <- merge_taxa_vec(
    ps, 
    group = clusters$cluster,
  )
  
  # save this out
  saveRDS(ps0, here(paste('data/sequencing_rpoB/phyloseq/myxococcus/clustered/ps_otu', percent_similarity[i], 'percent.rds', sep = '')))
  
  # do prevalence filtering
  # remove things that are only present in 3 or fewer samples & abundance > 100 overall
  ps0_sub <- microViz::tax_filter(ps0, min_prevalence = 4, min_total_abundance = 100)
  
  # save this out
  saveRDS(ps0_sub, here(paste('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu', percent_similarity[i], 'percent_filt.rds', sep = '')))
}