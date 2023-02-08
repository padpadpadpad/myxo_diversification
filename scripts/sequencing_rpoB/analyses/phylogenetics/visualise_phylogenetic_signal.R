# look at phylogenetic signal across levels of OTU similarity

#--------------------------#
# what this script does ####
#--------------------------#


# load in packages 
library(tidyverse)

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

d_dstat_sum <- group_by(d_dstat_boot, similarity2, habitat_preference) %>%
  summarise(mean = mean(estimated_d),
            .groups = 'drop')

# load in lambda
d_lambda <- map_df(files_lambda, read_csv, col_types = cols(similarity = col_character())) %>%
  mutate(similarity2 = ifelse(similarity == 'asv', 100, as.numeric(similarity)))
d_lambda_boot <- map_df(files_lambda_boot, read_csv, col_types = cols(similarity = col_character())) %>%
  mutate(similarity2 = ifelse(similarity == 'asv', 100, as.numeric(similarity)))

d_lambda_sum <- group_by(d_lambda_boot, similarity2, habitat_preference) %>%
  summarise(median = median(estimated_lambda),
            .groups = 'drop')

# make plot to see how D statistic changes across OTU cut-offs
ggplot() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(similarity2, estimated_d), position = position_jitter(width = 0.2), alpha = 0.2, d_dstat_boot) +
  geom_point(aes(similarity2, mean), size = 3, d_dstat_sum, shape =21, fill = 'white') +
  geom_line(aes(similarity2, estimated_d), d_dstat) +
  geom_point(aes(similarity2, estimated_d), size = 5, shape = 21, fill = 'white', d_dstat) +
  facet_wrap(~habitat_preference) +
  ylim(c(-0.7, 1)) +
  theme_bw() +
  labs(y = 'D statistic',
       x = 'OTU similarity')

# save out D statistic graph
ggsave('plots/sequencing_rpoB/analyses/d_statictic.png', height = 4, width = 10)

# make plot to see how lambda changes across OTU cut-offs
ggplot() +
  geom_point(aes(similarity2, estimated_lambda), position = position_jitter(width = 0.2), alpha = 0.2, d_lambda_boot) +
  geom_point(aes(similarity2, median), shape = 21, fill = 'white', size = 3, d_lambda_sum) +
  geom_line(aes(similarity2, estimated_lambda), d_lambda) +
  geom_point(aes(similarity2, estimated_lambda), size = 5, shape = 21, fill = 'white', d_lambda) +
  facet_wrap(~habitat_preference) +
  ylim(c(0, 1)) +
  theme_bw() +
  labs(y = "Pagel's lambda",
       x = 'OTU similarity')

# save out D statistic graph
ggsave('plots/sequencing_rpoB/analyses/lambda.png', height = 8, width = 10)

