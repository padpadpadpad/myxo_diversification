#---------------------#
# data exploration ####
#---------------------#

# looks at final 16S and myxo rpoB object and calculates summary stats

rm(list = ls())

# load packages ####
librarian::shelf(here, phyloseq, speedyseq, tidyverse)

here::i_am('scripts/misc/explore_16s_rpoB.R')

# set seed
set.seed(42)

# load in 16S data set
ps_16s <- readRDS(here('data/sequencing_16s/ps_16s_low_depth_removed.rds'))

# load in final myxo rpoB dataset
rpoB_seqs <- read_csv(here('data/sequencing_rpoB/phyloseq/myxococcus/myxo_seq_name_conversion.csv'))
summary(nchar(rpoB_seqs$seq))
ps_rpoB <- readRDS(here('data/sequencing_rpoB/phyloseq/ps_no_tree.rds')) %>%
  microViz::tax_filter(min_prevalence = 4, min_total_abundance = 100)

#----------------------#
# look at 16S reads ####
#----------------------#

# calculate number of reads that were myxococotta
d_16s <- psmelt(ps_16s) %>%
  janitor::clean_names() %>%
  mutate(myxo = ifelse(str_detect(phylum, 'Myxo'), 'yes', 'no'))

# check number of unique Myxo ASVs
filter(d_16s, myxo == 'yes') %>%
  pull(otu) %>%
  unique() %>%
  length()

# calculate summary stats for myxo reads per sample
d_16s_summary <- select(d_16s, sample, abundance, myxo) %>%
  #filter(abundance > 0) %>%
  group_by(sample) %>%
  mutate(total_16s = sum(abundance),
         present = ifelse(abundance > 0, 'yes', 'no')) %>%
  filter(myxo == 'yes') %>%
  group_by(sample, total_16s, present) %>%
  summarise(myxo_16s = sum(abundance),
            diversity_16s = n(),
            .groups = 'drop') %>%
  complete(sample, present, fill = list(myxo_16s = 0, diversity_16s = 0)) %>%
  group_by(sample) %>%
  mutate(total_16s = replace_na(total_16s, unique(na.omit(total_16s)))) %>%
  ungroup() %>%
  mutate(prop_16s = myxo_16s/total_16s) %>%
  filter(., present == 'yes')

#-----------------------#
# look at rpoB reads ####
#-----------------------#

# calculate number of reads that were myxococotta
d_rpoB <- psmelt(ps_rpoB) %>%
  janitor::clean_names() %>%
  mutate(myxo = ifelse(str_detect(phylum, 'Myxo'), 'yes', 'no'))

# check number of unique Myxo ASVs
filter(d_rpoB, myxo == 'yes') %>%
  pull(otu) %>%
  unique() %>%
  length()

# calculate summary stats for myxo reads per sample
d_rpoB_summary <- select(d_rpoB, sample, abundance, myxo) %>%
  #filter(abundance > 0) %>%
  group_by(sample) %>%
  mutate(total_rpoB = sum(abundance),
         present = ifelse(abundance > 0, 'yes', 'no')) %>%
  filter(myxo == 'yes') %>%
  group_by(sample, total_rpoB, present) %>%
  summarise(myxo_rpoB = sum(abundance),
            diversity_rpoB = n(),
            .groups = 'drop') %>%
  complete(sample, present, fill = list(myxo_rpoB = 0, diversity_rpoB = 0)) %>%
  group_by(sample) %>%
  mutate(total_rpoB = replace_na(total_rpoB, unique(na.omit(total_rpoB)))) %>%
  ungroup() %>%
  mutate(prop_rpoB = myxo_rpoB/total_rpoB) %>%
  filter(., present == 'yes')

left_join(d_16s_summary, d_rpoB_summary) %>%
  filter(!is.na(total_rpoB)) %>%
  mutate(fold_change = prop_rpoB/prop_16s) %>%
  pull(fold_change) %>%
  summary()
 