# visualise the habitat preference datasets

#--------------------------#
# what this script does ####
#--------------------------#

# read in each habitat preference dataset
# looks at how many ASVs are assigned to each habitat preference

# load packages
library(tidyverse)
library(MetBrewer)

# function to read in habitat preference dataset and put in percent similarity
read_habpref <- function(x){
  temp <- read.csv(x)
  # add column for filename
  name <- basename(x)
  temp$similarity = as.character(parse_number(name))
  if(str_detect(name, 'asv') == TRUE){temp$similarity <- 'asv'}
  
  return(temp)
}

# list all the habitat preference files
files <- list.files('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary', pattern = '.csv', full.names = TRUE)
files <- files[! basename(files) %in% c('habitat_preference_asv_new.csv', 'habitat_preference_asv_check.csv')]

# read them all in
d <- map_df(files, read_habpref)

# create a summary file for all of the OTU cut-offs
d_summary <- group_by(d, similarity, habitat_preference) %>%
  tally() %>%
  filter(similarity > 90) %>%
  group_by(similarity) %>%
  mutate(prop = n/sum(n),
         total_otu = sum(n),
         facet = paste(similarity, '% - total otus: ', total_otu, sep = ''),
         facet = ifelse(similarity == 'asv', paste('ASV - total otus: ', total_otu, sep = ''), facet),
         facet = forcats::fct_reorder(facet, similarity)) %>%
  ungroup() 

# set habitat colours
cols_hab <- met.brewer('Austria', n = 7)
names(cols_hab) <- c('marine_mud', 'freshwater', 'terrestrial', 'freshwater:terrestrial', 'freshwater:marine_mud:terrestrial', 'marine_mud:terrestrial', 'freshwater:marine_mud')
hab_labels <- c('marine specialist', 'freshwater specialist', 'land specialist', 'freshwater + land generalist', 'full generalist', 'marine + land generalist', 'freshwater + marine generalist')
names(hab_labels) <- c('marine_mud', 'freshwater', 'terrestrial', 'freshwater:terrestrial', 'freshwater:marine_mud:terrestrial', 'marine_mud:terrestrial', 'freshwater:marine_mud')

# create plot
ggplot(d_summary, aes(prop, forcats::fct_reorder(habitat_preference, n, .fun = mean))) +
  geom_col(aes(fill = habitat_preference)) +
  facet_wrap(~facet, labeller = labeller(facet = MicrobioUoE::letter_facets)) +
  theme_bw(base_size = 12) +
  geom_label(aes(label = n), size = MicrobioUoE::pts(10), nudge_x = 0.05, label.padding = unit(0.15, "lines")) +
  labs(x = 'Proportion of all species',
       y = 'Biome preference') +
  xlim(c(0,0.5)) +
  scale_fill_manual('Biome preference', values = cols_hab, labels = hab_labels) +
  scale_y_discrete(labels = hab_labels) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0))

# save this plot out
ggsave('plots/manuscript_plots/otu_similarity_diversity.png', last_plot(), height= 7, width =12)
