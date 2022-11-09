# look at proportion of transitions between states

# this script aims to look at how often transitions occur between and within habitat preferences

# load packages
library(tidyverse)
library(phytools)

# load in best output from simmap - the mac results because the model is better
simmap_res <- readRDS("data/sequencing_rpoB/processed/transition_rates/simmap_mac.rds")

# summarise number of switches between states and time spent in each state
simmap_summary <- describe.simmap(simmap_res, plot=FALSE)

# coerce transitions into dataframe
d_transitions <- as.data.frame(simmap_summary$count, col.names = colnames(simmap_summary$count)) %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(cols = c(everything(), -N, -iter), names_to = 'transition', values_to = 'n_transitions') %>%
  separate(transition, c('state_1', 'state_2'), sep = ',', remove = FALSE) %>%
  mutate(label = gsub(',', ' -> ', transition))

# calculate number of transitions from saline/non saline compared to within saline

d_marine <- mutate(d_transitions, marine_1 = ifelse(str_detect(state_1, 'mud_and_shore|generalist'), 'yes', 'no'),
                   marine_2 = ifelse(str_detect(state_2, 'mud_and_shore|generalist'), 'yes', 'no')) %>%
  group_by(., marine_1, marine_2, iter) %>%
  summarise(num_transitions = sum(n_transitions)) %>%
  group_by(iter) %>%
  mutate(prop = num_transitions/sum(num_transitions)) %>%
  ungroup() %>%
  group_by(marine_1, marine_2) %>%
  summarise(ave_prop = mean(prop),
            lower_ci = quantile(prop, 0.025),
            upper_ci = quantile(prop, 0.975),
            .groups = 'drop')

d_terrestrial <- mutate(d_transitions, land_1 = ifelse(str_detect(state_1, 'mud_and_shore|generalist'), 'yes', 'no'),
                   land_2 = ifelse(str_detect(state_2, 'terrestrial|generalist'), 'yes', 'no')) %>%
  group_by(., land_1, land_2, iter) %>%
  summarise(num_transitions = sum(n_transitions)) %>%
  group_by(iter) %>%
  mutate(prop = num_transitions/sum(num_transitions)) %>%
  ungroup() %>%
  group_by(land_1, land_2) %>%
  summarise(ave_prop = mean(prop),
            lower_ci = quantile(prop, 0.025),
            upper_ci = quantile(prop, 0.975),
            .groups = 'drop')

d_freshwater <- mutate(d_transitions, fresh_1 = ifelse(str_detect(state_1, 'mud_and_shore|generalist'), 'yes', 'no'),
                        fresh_2 = ifelse(str_detect(state_2, 'freshwater|generalist'), 'yes', 'no')) %>%
  group_by(., fresh_1, fresh_2, iter) %>%
  summarise(num_transitions = sum(n_transitions)) %>%
  group_by(iter) %>%
  mutate(prop = num_transitions/sum(num_transitions)) %>%
  ungroup() %>%
  group_by(fresh_1, fresh_2) %>%
  summarise(ave_prop = mean(prop),
            lower_ci = quantile(prop, 0.025),
            upper_ci = quantile(prop, 0.975),
            .groups = 'drop')

# So 76% of transitions are to and from habitat preferences that do not include marine mud, compared to 24% that do.
# We have quite a big split from marine mud to everything else. Consequently we can re-run our transition rate analysis but with saline/non-saline and fast / low diversification rate 
