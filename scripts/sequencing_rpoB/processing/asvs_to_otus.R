# Cluster existing phyloseq object to a specific OTU similarity.
# Also do prevalence filtering

#--------------------------#
# what this script does ####
#--------------------------#

# reads in myxo asv object
# does prevalence filtering: present in at least four samples and have a total abundance > 100
# aligns this object and creates a distance matrix
# clusters the ASVs at a given similarity
# saves out clustered phyloseq object

# load in packages
librarian::shelf(phyloseq, dada2, DECIPHER, here, tidyverse, mikemc/speedyseq)

# set number of processors to use
num_processors <- 4

# set where we are
here::i_am('scripts/sequencing_rpoB/processing/asvs_to_otus.R')

# read in phyloseq object
ps <- readRDS(here('data/sequencing_rpoB/phyloseq/myxococcus/clustered/ps_otu_asv.rds'))

# get taxonomy table from ps object
tax_table <- tax_table(ps) %>%
  data.frame() %>%
  rownames_to_column(var = 'tip_label') %>%
  janitor::clean_names()

# lets look at total abundance and total prevalence of each ASV
d_ps <- psmelt(ps) %>%
  janitor::clean_names() %>%
  mutate(presence = ifelse(abundance > 0, 1, 0)) %>%
  group_by(otu) %>%
  summarise(abundance = sum(abundance),
            prevalence = sum(presence),
            .groups = 'drop')

group_by(d_ps, prevalence) %>%
  tally()

ggplot(d_ps, aes(log10(abundance))) +
  geom_histogram()

# prevalence filter the raw ASV data
ps_sub <- microViz::tax_filter(ps, min_prevalence = 4, min_total_abundance = 100)

d_ps <- psmelt(ps_sub) %>%
  janitor::clean_names() %>%
  mutate(presence = ifelse(abundance > 0, 1, 0)) %>%
  group_by(otu) %>%
  summarise(abundance = sum(abundance),
            prevalence = sum(presence),
            .groups = 'drop')

group_by(d_ps, prevalence) %>%
  tally()

ggplot(d_ps, aes(log10(abundance))) +
  geom_histogram()

# save this out
saveRDS(ps_sub, here('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_asv_filt.rds'))

# read in myxo seq name conversion
seqs_df <- read.csv(here('data/sequencing_rpoB/phyloseq/myxococcus/myxo_seq_name_conversion.csv')) %>%
  filter(otu_name %in% taxa_names(ps_sub))

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
  
  # create new phyloseq object from the clusters
  ps0 <- merge_taxa_vec(
    ps_sub, 
    group = clusters$cluster,
  )

  # save this out
  saveRDS(ps0, here(paste('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_', percent_similarity[i], 'percent_filt.rds', sep = '')))
}