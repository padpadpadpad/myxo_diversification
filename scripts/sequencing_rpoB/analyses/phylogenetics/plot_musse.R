# results and plots of diversification analyses

#--------------------------#
# what this script does ####
#--------------------------#

# load in packages ####

librarian::shelf(ggtree, diversitree, tidybayes, ggdist, tidyverse)

# identify conflicts in the tidyverse packages and other packages
tidyverse_conflicts()

# read in datasets ####

# read in habitat preference
d_habpref <- read.csv(here('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_asv_new.csv'))

# read in phyloseq object and grab tax table
d_taxa <- readRDS(here('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_asv_filt.rds')) %>%
  phyloseq::tax_table() %>%
  data.frame() %>%
  janitor::clean_names() %>%
  rownames_to_column('otu')

# create d_meta
# use habitat_preference3 which collapses into marine mud generalists
d_meta <- left_join(select(d_habpref, otu, habitat_preference = habitat_preference3, num_present), select(d_taxa, otu:family))

# read in tree
tree <- read.tree(here('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_chronopl10.tre'))

# read in shift nodes
shiftnodes <- readRDS(here('data/sequencing_rpoB/processed/shiftnodes.rds'))

# read in colours for states
cols_hab <- readRDS(here('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_colours.rds'))

# setup for analyses ####

# reorder metadata to match tip labels of tree
d_meta <- tibble(tip_label = tree$tip.label) %>%
  left_join(., rename(d_meta, tip_label = otu))

# make sure order of habitat preference is the same in the trait vector as the tip labels
sum(d_meta$tip_label == tree$tip.label) == length(tree$tip.label)
# SUCCESS if TRUE

# use habitat preference 3 - which collapses marine mud generalists into a single state
unique(d_meta$habitat_preference)

hab_pref <- setNames(d_meta$habitat_preference, d_meta$tip_label)
hab_pref_num <- as.numeric(as.factor(hab_pref))
hab_pref_num <- setNames(hab_pref_num, d_meta$tip_label)

# coding from numeric to character
coding <- tibble(hab_pref = unname(hab_pref), hab_pref_num = unname(hab_pref_num)) %>%
  distinct() %>%
  arrange(hab_pref) %>%
  mutate(hab_pref_axis = gsub(':', '/ ', hab_pref),
         hab_pref_axis = gsub('_', ' ', hab_pref_axis),
         # rename the columns for easy naming
         initials = c('FG', 'FS', 'MG', 'MS', 'TS'))

coding

#--------------------------#
# musse model wrangling ####
#--------------------------#

# read in best markov model
fit_markov <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_mod_custom3_v2.rds')
max(fit_musse_no_se$par)*5

# read in musse models
fit_musse <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_musse.rds')
fit_musse_no_se <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_musse_no_se.rds')
fit_musse_no_sse <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_musse_no_sse.rds')
fit_musse_no_ss <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_musse_no_ss.rds')

# do model selection

# run anova
anova(no_sse = fit_musse_no_sse,
      no_ss = fit_musse_no_ss,
      no_se = fit_musse_no_se,
      fit_musse)

# make AIC table
d_musse_aic <- data.frame(model = c('full_sse', 'no_se', 'no_ss', 'no_sse'), aic = c(AIC(fit_musse), AIC(fit_musse_no_se), AIC(fit_musse_no_ss), AIC(fit_musse_no_sse))) %>%
  dplyr::arrange(., aic) %>%
  mutate(weights = round(MuMIn::Weights(aic), 3))

d_musse_aic
# model with no state specific extinction is favoured

# read in mcmc
musse_mcmc <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_mcmc_musse_no_se.rds')

# plot the speciation results
d_speciation <- select(musse_mcmc, i, starts_with('lambda')) %>%
  pivot_longer(cols = starts_with('lambda'), names_to = 'hab_pref_num', values_to = 'speciation_rate', names_prefix = 'lambda') %>%
  mutate(hab_pref_num = as.numeric(hab_pref_num)) %>%
  left_join(., coding)

d_speciation_summary <- group_by(d_speciation, hab_pref_num, hab_pref) %>%
  tidybayes::mean_qi(speciation_rate)

# grab out the maximum likelihood estimates
d_speciation_ml <- fit_musse_no_se$par[grepl('lambda', names(fit_musse_no_se$par))] %>%
  data.frame(speciation_rate = ., hab_pref_num = parse_number(names(.)), row.names = NULL) %>%
  left_join(., coding)

# plot speciation rate
ggplot(d_speciation_summary, aes(hab_pref, speciation_rate, col = hab_pref)) +
  stat_eye(aes(hab_pref, speciation_rate, col = hab_pref), d_speciation, show.legend = FALSE) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = cols_hab) +
  ylim(c(0,16)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  labs(x = 'Habitat preference',
       y = 'Speciation rate')

ggsave('plots/sequencing_rpoB/analyses/musse_speciation.png', last_plot(), height = 5, width = 6)

# how similar are the transition rates between the musse model and the simpler markov model
d_transitions_musse <- select(musse_mcmc, i, starts_with('q')) %>%
  pivot_longer(cols = starts_with('q'), names_to = 'transition', values_to = 'transition_rate', names_prefix = 'q') %>%
  mutate(state_1_num = as.numeric(substr(transition, 1,1)),
         state_2_num = as.numeric(substr(transition, 2,2))) %>%
  left_join(., select(coding, state_1_num = hab_pref_num, state_1 = hab_pref)) %>%
  left_join(., select(coding, state_2_num = hab_pref_num, state_2 = hab_pref)) %>%
  group_by(transition, state_1_num, state_1, state_2_num, state_2) %>%
  tidybayes::mean_qi(transition_rate)

d_transitions_markov <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_mcmc_custom3_v2.rds') %>%
  select(., i, starts_with('q')) %>%
  pivot_longer(cols = starts_with('q'), names_to = 'transition', values_to = 'transition_rate', names_prefix = 'q') %>%
  mutate(state_1_num = as.numeric(substr(transition, 1,1)),
         state_2_num = as.numeric(substr(transition, 2,2))) %>%
  left_join(., select(coding, state_1_num = hab_pref_num, state_1 = hab_pref)) %>%
  left_join(., select(coding, state_2_num = hab_pref_num, state_2 = hab_pref)) %>%
  group_by(transition, state_1_num, state_1, state_2_num, state_2) %>%
  tidybayes::mean_qi(transition_rate)

d_transitions <- left_join(
  select(d_transitions_musse, transition, state_1_num, state_1, state_2_num, state_2, rate_musse = transition_rate, lower_musse = .lower, upper_musse = .upper),
  select(d_transitions_markov, transition, state_1_num, state_1, state_2_num, state_2, rate_markov = transition_rate, lower_markov = .lower, upper_markov = .upper)
)

ggplot(d_transitions, aes(rate_markov, rate_musse)) +
  geom_point() +
  geom_linerange(aes(ymin = lower_musse, ymax = upper_musse)) +
  geom_linerange(aes(xmin = lower_markov, xmax = upper_markov)) +
  xlim(c(0, 60)) +
  ylim(c(0, 65))

# look at which rates correlate with each other from the Musse model
select(musse_mcmc, i, starts_with('q')) %>%
  slice_sample(n = 250) %>%
  ggpairs(., columns = 2:ncol(.),
          lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw()
ggsave('plots/sequencing_rpoB/analyses/discrete_character_evolution/pairs_plot_musse.png', last_plot(), height = 12, width = 14)
