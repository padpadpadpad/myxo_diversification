# look at phylogenetic signal across levels of OTU similarity

#--------------------------#
# what this script does ####
#--------------------------#

# 1. reads in estimates of phylogenetic signal at different otu cut-offs
# 2. creates plot to visualise these

# load in packages 
librarian::shelf(tidyverse, patchwork)

# list files

# d statistic
files_d <- list.files('data/sequencing_rpoB/processed/phylogenetic_signal', full.names = TRUE, pattern = 'd_statistic')
files_d_boot <- files_d[str_detect(files_d, 'boot')]
files_d <- files_d[!files_d %in% files_d_boot]

# lambda
files_lambda <- list.files('data/sequencing_rpoB/processed/phylogenetic_signal', full.names = TRUE, pattern = 'lambda')
files_lambda_boot <- files_lambda[str_detect(files_lambda, 'boot')]
files_lambda <- files_lambda[!files_lambda %in% files_lambda_boot]

# load in D statistic
d_dstat <- map_df(files_d, read_csv, col_types = cols(similarity = col_character())) %>%
  mutate(similarity2 = ifelse(similarity == 'asv', 100, as.numeric(similarity)))
d_dstat_boot <- map_df(files_d_boot, read_csv, col_types = cols(similarity = col_character())) %>%
  mutate(similarity2 = ifelse(similarity == 'asv', 100, as.numeric(similarity)))

# summaries the bootstrap estimates
d_dstat_sum <- group_by(d_dstat_boot, similarity2, habitat_preference) %>%
  summarise(mean = mean(estimated_d),
            .groups = 'drop')

# load in lambda
d_lambda <- map_df(files_lambda, read_csv, col_types = cols(similarity = col_character())) %>%
  mutate(similarity2 = ifelse(similarity == 'asv', 100, as.numeric(similarity)))
d_lambda_boot <- map_df(files_lambda_boot, read_csv, col_types = cols(similarity = col_character())) %>%
  mutate(similarity2 = ifelse(similarity == 'asv', 100, as.numeric(similarity)))

# summarise the lambda bootstraps
d_lambda_sum <- group_by(d_lambda_boot, similarity2, habitat_preference) %>%
  summarise(median = median(estimated_lambda),
            .groups = 'drop')

# make plot to see how D statistic changes across OTU cut-offs
p_1 <- ggplot() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point(aes(similarity2, estimated_d, col = habitat_preference), position = position_jitter(width = 0.2), alpha = 0.2, d_dstat_boot, show.legend = FALSE) +
  geom_point(aes(similarity2, mean), size = 2, d_dstat_sum, shape =21, fill = 'white') +
  geom_line(aes(similarity2, estimated_d), d_dstat) +
  geom_point(aes(similarity2, estimated_d), size = 3, shape = 21, fill = 'white', d_dstat) +
  facet_wrap(~habitat_preference, labeller = labeller(.default = MicrobioUoE::letter_facets, habitat_preference = c(freshwater = 'freshwater (yes or no)', marine_mud = 'marine (yes or no)', terrestrial = 'land (yes or no)'))) +
  ylim(c(-0.7, 1)) +
  theme_bw(base_size = 14) +
  labs(y = 'D statistic',
       x = 'OTU similarity (%)') +
  scale_color_manual('Habitat cluster', values = c('#5ECAE2', '#3911EE', '#53C20A'), labels = c('Freshwater (yes or no)', 'Marine', 'Land (yes or no)')) +
  scale_x_continuous(breaks = c(91:100), minor_breaks = NULL) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10))

# make plot to see how lambda changes across OTU cut-offs
p_2 <- ggplot() +
  geom_point(aes(similarity2, estimated_lambda), position = position_jitter(width = 0.2), alpha = 0.2, filter(d_lambda_boot, habitat_preference == 'all')) +
  geom_point(aes(similarity2, median), shape = 21, fill = 'white', size = 2, filter(d_lambda_sum, habitat_preference == 'all')) +
  geom_line(aes(similarity2, estimated_lambda), filter(d_lambda, habitat_preference == 'all')) +
  geom_point(aes(similarity2, estimated_lambda), size = 5, shape = 21, fill = 'white', filter(d_lambda, habitat_preference == 'all')) +
  ylim(c(0, 1)) +
  theme_bw(base_size = 14) +
  labs(y = "Pagel's lambda",
       x = 'OTU similarity (%)',
       title = '(d) 7-state biome preference') +
  scale_x_continuous(breaks = c(91:100), minor_breaks = NULL) +
  theme(plot.title = element_text(size = 12))

p_1 / p_2 + plot_layout(heights = c(0.4, 0.6))

# save out D statistic graph
ggsave('plots/manuscript_plots/phylogenetic_signal.png', height = 7, width = 8)

