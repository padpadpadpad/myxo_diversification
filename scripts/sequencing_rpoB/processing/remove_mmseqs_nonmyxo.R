# Script to remove any ambiguous taxa from the raw Myxococcus phyloseq object

#--------------------------#
# what this script does ####
#--------------------------#

# reads in the Myxococcota phyloseq object
# reads in the output from MMSeqs LCA algorithm
# removes ASVs which are not assigned to Myxococcus using the LCA algorithm
# looks at similarity between MMSeqs and AssignTaxonomy approach
# saves out new Myxococcota phyloseq obect which is then used

# load packages
library(phyloseq)
library(here)
library(tidyverse)

# set where we are
here::i_am('scripts/sequencing_rpoB/processing/remove_mmseqs_nonmyxo.R')

# read in phyloseq object
ps <- readRDS(here('data/sequencing_rpoB/phyloseq/myxococcus/ps_phyloseq_myxo.rds'))

ps 

# grab taxonomy table from phyloseq object
d_taxa <- tax_table(ps) %>%
  data.frame() %>%
  rownames_to_column('otu')

# read in MMSeqs2 taxonomy dataset

# we ran these commands to come up with the taxonomic assignments using MMSeqs2
# mmseqs taxonomy seqTaxDB taxonomy/gtdb_r202_rpoB_mmseqs taxonomyResult_lca3_95 tmpFolder --search-type 3 --tax-lineage 1 --majority 0.95 --vote-mode 1 --lca-mode 3 --orf-filter 0
# mmseqs createtsv seqTaxDB taxonomyResult_lca3_95 taxonomy_assignments_lca3_95.tsv

d_taxa_mmseq <- read_tsv('data/sequencing_rpoB/phyloseq/myxococcus/mmseqs_lca3_95.tsv', col_names = FALSE)
names(d_taxa_mmseq) <- c('otu', 'id', 'taxon_rank', 'name', 'taxonomy')

# wrangle MMSeqs object
d_taxa_mmseq <- mutate(d_taxa_mmseq, taxonomy = gsub('d_|p_|c_|o_|f_|g_|s_', '', taxonomy)) %>%
  separate(., taxonomy, c('Kingdom', 'Phylum', 'Class', 'Order','Family','Genus', 'Species'), sep = ';') %>%
  select(otu, Kingdom:Species)

# how many of the OTUs are not assigned to myxococotta using MMSeqs
to_remove <- filter(d_taxa_mmseq, !str_detect(Phylum, 'Myxococcota') | is.na(Phylum))
# 76 ASVs

# remove these from the dataset
to_keep <- taxa_names(ps)[!taxa_names(ps) %in% to_remove$otu]
ps_sub <- prune_taxa(to_keep, ps)

# compare similarity between levels of taxonomy between MMSeqs and assignTaxonomy

# filter mmseqs and phyloseq taxonomies to just keep ones which are all myxo
d_taxa <- filter(d_taxa, otu %in% to_keep)
d_taxa_mmseq <- filter(d_taxa_mmseq, otu %in% to_keep)

# how similar are the taxonomies
d_compare <- left_join(d_taxa %>% pivot_longer(Kingdom:Species, names_to = 'rank', values_to = 'dada2'), 
                       pivot_longer(d_taxa_mmseq, Kingdom:Species, names_to = 'rank', values_to = 'mmseq')) %>%
  group_by(rank) %>%
  summarise(total = n(),
            na_dada2 = sum(is.na(dada2)),
            match = sum(dada2 == mmseq, na.rm = TRUE),
            mismatch = sum(dada2 != mmseq, na.rm = TRUE),
            na_mmseq = sum(is.na(mmseq)),
            na_mmseq_unique = sum(is.na(mmseq) & !is.na(dada2)),
            prop_match_mismatch = match/(match + mismatch)) %>%
  arrange(na_dada2)

# generally MMSeqs2 matches up really well with assignTaxonomy
# proportion of mismatches is really low

# so we used MMSeqs2 to assign taxonomy from now on
# replace the tax table of the phyloseq and resave it

# save out new phyloseq object

# convert MMSeqs2 object to OTU table for phyloseq
new_tax_table <- column_to_rownames(d_taxa_mmseq, var = 'otu') %>%
  as.matrix() %>%
  tax_table()

# replace tax table with new one from MMSeqs2
tax_table(ps_sub) <- new_tax_table

saveRDS(ps_sub, here('data/sequencing_rpoB/phyloseq/myxococcus/clustered/ps_otu_asv.rds'))
