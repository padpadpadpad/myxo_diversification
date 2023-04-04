#---------------------------------------------------------#
# script to run models of discrete character evolution ####
#---------------------------------------------------------#

# what this script does ####

# 1.
# 2.

# load packages ####
librarian::shelf(here, tidyverse, ggtree, ggnewscale, RColorBrewer, patchwork, phytools, ggpp, castor, diversitree, tidygraph, igraph, ggraph, GGally, ggrepel, flextable, ggridges, hisse, padpadpadpad/MicrobioUoE)

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

#custom_function ####

# function for getting a data frame from a diversitree object
get_diversitree_df <- function(div_obj, trait_vec, replace_vec){
  
  if(is.null(div_obj$par.full)){div_obj$par.full <- div_obj$par}
  
  temp <- tibble(param = names(div_obj$par.full)) %>%
    mutate(state_1_num = as.numeric(substr(param, 2,2)),
         state_2_num = as.numeric(substr(param, 3,3)),
         transition_rate = unlist(div_obj$par.full),
         state_1 = stringi::stri_replace_all_regex(state_1_num, pattern = trait_vec, replacement = replace_vec, vectorize=FALSE),
         state_2 = stringi::stri_replace_all_regex(state_2_num, pattern = trait_vec, replacement = replace_vec, vectorize=FALSE)) %>%
    select(param, state_1, state_2, state_1_num, state_2_num, transition_rate) %>%
    mutate(free_param = ifelse(param %in% names(div_obj$par), 'yes', 'no'),
           num_params = length(div_obj$par))
  
  return(temp)
  
}

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

#-------------------------------------------------------------#
# run diversitree analysis of discrete character evolution ####
#-------------------------------------------------------------#

# all of these methods needs a likelihood function, we can build a Mkn model
lik_ard <- make.mkn(tree, hab_pref_num, k = max(hab_pref_num))

argnames(lik_ard)

# make symmetric rates model
lik_sym <- constrain(lik_ard, 
                     q12~q21, q13~q31, q14~q41, q15~q51,
                     q23~q32, q24~q42, q25~q52,
                     q34~q43, q35~q53,
                     q45~q54)

argnames(lik_sym)

# make equal rates model
lik_er <- constrain(lik_ard,
                    q13~q12, q14~q12, q15~q12,
                    q21~q12, q23~q12, q24~q12, q25~q12,
                    q31~q12, q32~q12, q34~q12, q35~q12,
                    q41~q12, q42~q12, q43~q12, q45~q12,
                    q51~q12, q52~q12, q53~q12, q54~q12)

argnames(lik_er)

lik_transient <- constrain(lik_ard,
                           q24~0, q25~0, q42~0, q52~0, q14~0, q41~0, q45~0, q54~0)

# need to pass start values to it - can grab these from the ape::ace, but we will just pass an average rate to the model
inits_ard <- rep(1, length(argnames(lik_ard)))
inits_sym <- rep(1, length(argnames(lik_sym)))
inits_er <- rep(1, length(argnames(lik_er)))
inits_trans <- rep(1, length(argnames(lik_transient)))

# run first three models
mod_er <- find.mle(lik_er, inits_er, method = 'subplex', control = list(maxit = 50000))
mod_ard <- find.mle(lik_ard, inits_ard, method = 'subplex', control = list(maxit = 50000))
mod_trans <- find.mle(lik_transient, inits_trans, method = 'subplex', control = list(maxit = 50000))
mod_sym <- find.mle(lik_sym, inits_sym, method = 'subplex', control = list(maxit = 500))

# compare models
AIC(mod_sym, mod_ard, mod_er, mod_trans) %>% arrange(AIC)

# sort parameter estimates
mod_ard$par %>% sort()

# get parameters close to 0.
mod_ard$par[mod_ard$par < 1e-03] %>% names(.)

# make custom matrix model
lik_custom1 <- constrain(lik_ard, 
                     q14~0, q24~0, q25~0, q41~0, q53~0)

# make start parameters
inits_custom1 <- rep(1, length(argnames(lik_custom1)))

# refit model
mod_custom1 <- find.mle(lik_custom1, inits_custom1, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1) %>% arrange(AIC)

# filter for only estimated parameters
mod_custom1$par.full[names(mod_custom1$par.full) %in% argnames(lik_custom1)] 

# sort parameter estimates
mod_custom1$par %>% sort()

# make custom matrix model again
lik_custom2 <- constrain(lik_ard, 
                         q14~0, q24~0, q25~0, q41~0, q53~0,
                         q42~0)

# make start parameters
inits_custom2 <- rep(1, length(argnames(lik_custom2)))

# run model
mod_custom2 <- find.mle(lik_custom2, inits_custom2, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2) %>% arrange(AIC)

# sort parameter estimates
mod_custom2$par %>% sort()

# make custom matrix model
lik_custom3 <- constrain(lik_ard, 
                         q14~0, q24~0, q25~0, q41~0, q53~0,
                         q42~0,
                         q45~0)

# make start parameters
inits_custom3 <- rep(1, length(argnames(lik_custom3)))

# run model
mod_custom3 <- find.mle(lik_custom3, inits_custom3, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3) %>% arrange(AIC)

# sort parameter estimates
mod_custom3$par %>% sort()

# make custom matrix model
lik_custom4 <- constrain(lik_ard, 
                         q14~0, q24~0, q25~0, q41~0, q53~0,
                         q42~0,
                         q45~0,
                         q54~0)

# make start parameters
inits_custom4 <- rep(1, length(argnames(lik_custom4)))

# run model
mod_custom4 <- find.mle(lik_custom4, inits_custom4, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3, mod_custom4) %>% arrange(AIC)

# sort parameter estimates
mod_custom4$par %>% sort()

# make custom matrix model
lik_custom5 <- constrain(lik_ard, 
                         q14~0, q24~0, q25~0, q41~0, q53~0,
                         q42~0,
                         q45~0,
                         q54~0,
                         q52~0)

# make start parameters
inits_custom5 <- rep(1, length(argnames(lik_custom5)))

# run model
mod_custom5 <- find.mle(lik_custom5, inits_custom5, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3, mod_custom4, mod_custom5) %>% arrange(AIC)

# sort parameter estimates
mod_custom5$par %>% sort()

# make custom matrix model
lik_custom6 <- constrain(lik_ard, 
                         q14~0, q24~0, q25~0, q41~0, q53~0,
                         q42~0,
                         q45~0,
                         q54~0,
                         q52~0,
                         q43~0)

# make start parameters
inits_custom6 <- rep(1, length(argnames(lik_custom6)))

# run model
mod_custom6 <- find.mle(lik_custom6, inits_custom6, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3, mod_custom4, mod_custom5, mod_custom6) %>% arrange(AIC)

# lets look at model weights
d_aic <- AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3, mod_custom4, mod_custom5, mod_custom6) %>%
  data.frame() %>%
  mutate(log_lik = c(mod_sym$lnLik, mod_ard$lnLik, mod_er$lnLik, mod_trans$lnLik, mod_custom1$lnLik, mod_custom2$lnLik, mod_custom3$lnLik, mod_custom4$lnLik, mod_custom5$lnLik, mod_custom6$lnLik),
         ntips = 2621) %>%
  janitor::clean_names() %>%
  rownames_to_column(var = 'model') %>%
  mutate(aicc = -2*log_lik + 2*df*(ntips/(ntips - df - 1)),
         aic_weight = round(MuMIn::Weights(aic), 2),
         aicc_weight = round(MuMIn::Weights(aicc), 2)) %>%
  arrange(-aic_weight)
# mod custom 5 is favoured pretty highly

# save out mcmc chains
saveRDS(mod_custom5, 'data/sequencing_rpoB/processed/transition_rates/mod_custom_5.rds')

# plot_transition_matrix
diversitree_df <- get_diversitree_df(mod_custom5, coding$hab_pref_num, coding$hab_pref)

diversitree_df %>%
  left_join(., select(coding, state_1 = hab_pref, state_1_num = hab_pref_num, state_1_label = hab_pref_axis)) %>%
  left_join(., select(coding, state_2 = hab_pref, state_2_num = hab_pref_num, state_2_label = hab_pref_axis)) %>%
  mutate(transition_rate = round(transition_rate, 2)) %>%
  ggplot(., aes(forcats::fct_reorder(state_2_label, state_2_num), forcats::fct_reorder(state_1_label, desc(state_1_num)))) +
  geom_tile(aes(alpha = transition_rate, col = free_param), width = 0.9, height = 0.9, size = 1.1) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
  legend.position = 'none',
  axis.text.x.top = element_text(angle = 90, vjust = 0.5),
  plot.title.position = "plot") +
  scale_alpha_continuous(range = c(0, 0.6)) +
  geom_text(aes(label = transition_rate), size = MicrobioUoE::pts(10)) +
  scale_x_discrete(position = 'top', labels = scales::label_wrap(13)) +
  scale_y_discrete(position = 'left', labels = scales::label_wrap(13)) +
  labs(y = 'From',
  x = 'To',
  title = paste('all rates different with', length(mod_custom5$par), 'free parameters', sep = ' ')) +
  coord_fixed() +
  scale_color_manual(values = c('red', 'black'))

ggsave('plots/sequencing_rpoB/analyses/discrete_character_evolution/transition_matrix.png', last_plot(), height = 5, width = 7)

#-----------------------------------------------------------------#
# look at uncertainty in these model estimates by running MCMC ####
#-----------------------------------------------------------------#

# set up initial start values
inits_mcmc <- mod_custom5$par

## # set up upper and lower limits - limit the values to be <3 times the max value
lower_mcmc <- rep(0, length(inits_mcmc))
upper_mcmc <- rep(max(inits_mcmc)*10, length(inits_mcmc))

# run first mcmc to tune w
fit_mcmc <- mcmc(lik_custom5, inits_mcmc, nsteps = 10, w = 0.1, upper = upper_mcmc, lower = lower_mcmc)

# tune w for each parameter
w <- diff(sapply(fit_mcmc[2:(ncol(fit_mcmc)-1)], quantile, c(.05, .95)))

# run second mcmc to tune w
fit_mcmc2 <- mcmc(lik_custom5, inits_mcmc, nsteps=100, w=w, upper = upper_mcmc, lower = lower_mcmc)

## # tune w for each parameter
w <- diff(sapply(fit_mcmc2[2:(ncol(fit_mcmc2)-1)], quantile, c(.05, .95)))

# run third mcmc for 1000 iter
fit_mcmc3 <- mcmc(lik_custom5, inits_mcmc, nsteps=1000, w=w, upper = upper_mcmc, lower = lower_mcmc)

# save out mcmc chains
saveRDS(fit_mcmc3, 'data/sequencing_rpoB/processed/transition_rates/asv_mcmc_trans.rds')

# make data long format
d_mcmc <- pivot_longer(fit_mcmc3, names_to = 'param', values_to = 'transition_rate', cols = starts_with('q')) %>%
  left_join(., select(diversitree_df, param, state_1, state_2)) %>%
  left_join(., select(coding, state_1 = hab_pref, state_1_num = hab_pref_num, state_1_label = initials)) %>%
  left_join(., select(coding, state_2 = hab_pref, state_2_num = hab_pref_num, state_2_label = initials)) %>%
  mutate(parameter = paste(state_1_label, '->', state_2_label))

# find 95% CIs and bind with ML estimates
d_mcmc_summary <- d_mcmc %>%
  group_by(parameter, param, state_1, state_2) %>%
  tidybayes::mean_qi(transition_rate) %>%
  left_join(., select(diversitree_df, state_1, state_2, ml_estimate = transition_rate))

# plot ridge plot
ggplot(d_mcmc, aes(transition_rate, forcats::fct_reorder(parameter, transition_rate))) +
  geom_density_ridges() +
  theme_bw(base_size = 14) +
  labs(x = 'transition rate',
       y = 'transition')

ggplot(d_mcmc_summary, aes(transition_rate, forcats::fct_reorder(parameter, transition_rate))) +
  geom_linerange(aes(xmin = .lower, xmax = .upper)) +
  geom_point(size = 3) +
  geom_point(aes(x = ml_estimate), size = 3, col = 'red') +
  theme_bw(base_size = 14) +
  labs(x = 'transition rate',
       y = 'transition',
       caption = 'red points are ML estimate\nblack points are MCMC average')

# look at which rates correlate with each
d_mcmc %>%
  select(parameter, transition_rate, i) %>%
  pivot_wider(names_from = parameter, values_from = transition_rate) %>%
  slice_sample(n = 250) %>%
  ggpairs(., columns = 2:ncol(.),
          lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw()

ggsave('plots/sequencing_rpoB/analyses/discrete_character_evolution/pairs_plot.png', last_plot(), height = 12, width = 14)

#------------------------------------#
# plot best transition rate model ####
#------------------------------------#

label_wrap2 <- function(x, width){
  unlist(lapply(strwrap(x, width = width, simplify = FALSE), 
                paste0, collapse = "\n"))
}

# calculate proportion of number of species each state has
d_habpref_summary <- group_by(d_meta, habitat_preference) %>%
  tally() %>%
  ungroup() %>%
  mutate(prop = n/sum(n))

# set colours
cols_hab <- readRDS(here('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_colours.rds'))

# make very simple plot to grab legend from
p_legend <- ggplot(d_habpref_summary, aes(habitat_preference, prop, col = habitat_preference)) +
  geom_point() +
  scale_color_manual('Habitat preference', values = cols_hab) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size=5)))
p_legend <- cowplot::get_legend(p_legend)

# turn transition matrix into network to plot
d_network <- as_tbl_graph(select(diversitree_df, state_1, state_2, transition_rate)) %>%
  activate(edges) %>%
  filter(!is.na(transition_rate) & transition_rate > 0) %>%
  activate(nodes) %>%
  left_join(., select(d_habpref_summary, name = habitat_preference, prop, n)) %>%
  left_join(., tibble(name = names(cols_hab))) %>%
  mutate(order = c(2, 1, 5, 4, 3)) %>%
  arrange(order)

p <- ggraph(d_network, layout = 'linear', circular = TRUE) + 
  geom_edge_fan(aes(alpha = transition_rate, 
                    width = transition_rate,
                    label = round(transition_rate, 2)),
                arrow = arrow(length = unit(4, 'mm')),
                end_cap = circle(10, 'mm'),
                start_cap = circle(10, 'mm'),
                angle_calc = 'along',
                label_dodge = unit(2.5, 'mm'),
                strength = 1.3) + 
  geom_node_point(aes(size = prop,
                      col = name)) +
  theme_void() +
  #geom_node_label(aes(label = label, x=xmin), repel = TRUE) +
  scale_size(range = c(5,20)) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_color_manual('Habitat preference', values = cols_hab) +
  scale_fill_manual('Habitat preference', values = cols_hab)

# grab data for points
point_data <- p$data %>%
  select(x, y, name) %>%
  mutate(nudge_x = ifelse(x < 0, -0.2, 0.2),
         nudge_y = ifelse(y < 0, -0.2, 0.2))

p_best <- p + 
  #geom_label(aes(nudge_x + x, nudge_y+y, label = label_wrap2(name, 15)), point_data, size = MicrobioUoE::pts(12)) +
  #coord_cartesian(clip = "off") +
  xlim(c(min(point_data$x) + min(point_data$nudge_x)), max(point_data$x) + max(point_data$nudge_x)) +
  ylim(c(min(point_data$y) + min(point_data$nudge_y)), max(point_data$y) + max(point_data$nudge_y)) +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = 'white', colour = 'white'))

p_best + (wrap_elements(p_legend)/plot_spacer()) + plot_layout(widths = c(0.8, 0.2))

# save out model
ggsave(here('plots/sequencing_rpoB/analyses/discrete_character_evolution/transition_plot_diversitree.png'), last_plot(), height = 5.5, width = 7)

# plot other models based on AIC weights

# turn transition matrix into network to plot
d_network <- get_diversitree_df(mod_custom3, coding$hab_pref_num, coding$hab_pref) %>%
  select(., state_1, state_2, transition_rate) %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  filter(!is.na(transition_rate) & transition_rate > 0) %>%
  activate(nodes) %>%
  left_join(., select(d_habpref_summary, name = habitat_preference, prop, n)) %>%
  left_join(., tibble(name = names(cols_hab))) %>%
  mutate(order = c(2, 1, 5, 4, 3)) %>%
  arrange(order)

p <- ggraph(d_network, layout = 'linear', circular = TRUE) + 
  geom_edge_fan(aes(alpha = transition_rate, 
                    width = transition_rate,
                    label = round(transition_rate, 2)),
                arrow = arrow(length = unit(4, 'mm')),
                end_cap = circle(10, 'mm'),
                start_cap = circle(10, 'mm'),
                angle_calc = 'along',
                label_dodge = unit(2.5, 'mm'),
                strength = 1.3) + 
  geom_node_point(aes(size = prop,
                      col = name)) +
  theme_void() +
  #geom_node_label(aes(label = label, x=xmin), repel = TRUE) +
  scale_size(range = c(5,20)) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_color_manual('Habitat preference', values = cols_hab) +
  scale_fill_manual('Habitat preference', values = cols_hab)

# grab data for points
point_data <- p$data %>%
  select(x, y, name) %>%
  mutate(nudge_x = ifelse(x < 0, -0.2, 0.2),
         nudge_y = ifelse(y < 0, -0.2, 0.2))

p_2 <- p + 
  #geom_label(aes(nudge_x + x, nudge_y+y, label = label_wrap2(name, 15)), point_data, size = MicrobioUoE::pts(12)) +
  #coord_cartesian(clip = "off") +
  xlim(c(min(point_data$x) + min(point_data$nudge_x)), max(point_data$x) + max(point_data$nudge_x)) +
  ylim(c(min(point_data$y) + min(point_data$nudge_y)), max(point_data$y) + max(point_data$nudge_y)) +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = 'white', colour = 'white')) +
  labs(title = '(b) Simplified ARD model 3',
       subtitle = 'AIC model weight = 0.19')

# turn transition matrix into network to plot
d_network <- get_diversitree_df(mod_custom4, coding$hab_pref_num, coding$hab_pref) %>%
  select(., state_1, state_2, transition_rate) %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  filter(!is.na(transition_rate) & transition_rate > 0) %>%
  activate(nodes) %>%
  left_join(., select(d_habpref_summary, name = habitat_preference, prop, n)) %>%
  left_join(., tibble(name = names(cols_hab))) %>%
  mutate(order = c(2, 1, 5, 4, 3)) %>%
  arrange(order)

p <- ggraph(d_network, layout = 'linear', circular = TRUE) + 
  geom_edge_fan(aes(alpha = transition_rate, 
                    width = transition_rate,
                    label = round(transition_rate, 2)),
                arrow = arrow(length = unit(4, 'mm')),
                end_cap = circle(10, 'mm'),
                start_cap = circle(10, 'mm'),
                angle_calc = 'along',
                label_dodge = unit(2.5, 'mm'),
                strength = 1.3) + 
  geom_node_point(aes(size = prop,
                      col = name)) +
  theme_void() +
  #geom_node_label(aes(label = label, x=xmin), repel = TRUE) +
  scale_size(range = c(5,20)) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_color_manual('Habitat preference', values = cols_hab) +
  scale_fill_manual('Habitat preference', values = cols_hab)

# grab data for points
point_data <- p$data %>%
  select(x, y, name) %>%
  mutate(nudge_x = ifelse(x < 0, -0.2, 0.2),
         nudge_y = ifelse(y < 0, -0.2, 0.2))

p_3 <- p + 
  #geom_label(aes(nudge_x + x, nudge_y+y, label = label_wrap2(name, 15)), point_data, size = MicrobioUoE::pts(12)) +
  #coord_cartesian(clip = "off") +
  xlim(c(min(point_data$x) + min(point_data$nudge_x)), max(point_data$x) + max(point_data$nudge_x)) +
  ylim(c(min(point_data$y) + min(point_data$nudge_y)), max(point_data$y) + max(point_data$nudge_y)) +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = 'white', colour = 'white')) +
  labs(title = '(c) Simplified ARD model 4',
       subtitle = 'AIC model weight = 0.15')

# turn transition matrix into network to plot
d_network <- get_diversitree_df(mod_trans, coding$hab_pref_num, coding$hab_pref) %>%
  select(., state_1, state_2, transition_rate) %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  filter(!is.na(transition_rate) & transition_rate > 0) %>%
  activate(nodes) %>%
  left_join(., select(d_habpref_summary, name = habitat_preference, prop, n)) %>%
  left_join(., tibble(name = names(cols_hab))) %>%
  mutate(order = c(2, 1, 5, 4, 3)) %>%
  arrange(order)

p <- ggraph(d_network, layout = 'linear', circular = TRUE) + 
  geom_edge_fan(aes(alpha = transition_rate, 
                    width = transition_rate,
                    label = round(transition_rate, 2)),
                arrow = arrow(length = unit(4, 'mm')),
                end_cap = circle(10, 'mm'),
                start_cap = circle(10, 'mm'),
                angle_calc = 'along',
                label_dodge = unit(2.5, 'mm'),
                strength = 1.3) + 
  geom_node_point(aes(size = prop,
                      col = name)) +
  theme_void() +
  #geom_node_label(aes(label = label, x=xmin), repel = TRUE) +
  scale_size(range = c(5,20)) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_color_manual('Habitat preference', values = cols_hab) +
  scale_fill_manual('Habitat preference', values = cols_hab)

# grab data for points
point_data <- p$data %>%
  select(x, y, name) %>%
  mutate(nudge_x = ifelse(x < 0, -0.2, 0.2),
         nudge_y = ifelse(y < 0, -0.2, 0.2))

p_4 <- p + 
  #geom_label(aes(nudge_x + x, nudge_y+y, label = label_wrap2(name, 15)), point_data, size = MicrobioUoE::pts(12)) +
  #coord_cartesian(clip = "off") +
  xlim(c(min(point_data$x) + min(point_data$nudge_x)), max(point_data$x) + max(point_data$nudge_x)) +
  ylim(c(min(point_data$y) + min(point_data$nudge_y)), max(point_data$y) + max(point_data$nudge_y)) +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = 'white', colour = 'white')) +
  labs(title = '(d) Transient model',
       subtitle = 'AIC model weight = 0.14')

# make layout
design <- "
  123
  453
"

p_best <- p_best + labs(title = '(a) Simplified ARD model 5',
                        subtitle = 'AIC model weight = 0.4')

p_best + p_2 + p_legend + p_3 + p_4 + plot_layout(design = design, widths = c(0.4, 0.4, 0.2))

# save out model
ggsave(here('plots/sequencing_rpoB/analyses/discrete_character_evolution/transition_plot_multiple.png'), last_plot(), height = 9, width = 12)

#------------------------------------------------#
# look at whether states are a source or sink ####
#------------------------------------------------#

d_source_sink_rate <- select(diversitree_df, away = state_1, into = state_2, transition_rate) %>%
  pivot_longer(cols = c(away, into), names_to = 'direction', values_to = 'habitat_preference') %>%
  group_by(habitat_preference, direction) %>%
  summarise(total_rate = sum(transition_rate), .groups = 'drop') %>%
  pivot_wider(names_from = direction, values_from = total_rate) %>%
  mutate(source_sink1 = away / into,
         source_sink2 = into - away)

table_rate <- select(d_source_sink_rate, habitat_preference, away, into, source_sink1) %>%
  mutate(across(away:source_sink1, ~round(.x, 2)),
         habitat_preference = gsub(':', ' + ', habitat_preference),
         habitat_preference = gsub('_', ' ', habitat_preference)) %>%
  arrange(desc(source_sink1)) %>%
  flextable(.) %>%
  set_header_labels(habitat_preference = 'habitat preference',
                    source_sink1 = 'source sink ratio') %>%
  align(align = 'center', part = 'all') %>%
  align(align = 'left', part = 'body', j = 1) %>%
  bold(part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit()

# save out
save_as_image(table_rate, here('plots/sequencing_rpoB/analyses/discrete_character_evolution/source_sink_rate.png'), zoom = 3, webshot = 'webshot2')

#-------------------------#
# CODE NOT RAN ANYMORE ####
#-------------------------#

#-------------------------------------#
# run stochastic character mapping ####
#-------------------------------------#

# create transition matrix
num_states <- unique(hab_pref) %>% length()

# we will set up the custom matrix, this has 49 numbers and then we have to set the correct numbers to 0
best_matrix <- matrix(1:num_states^2, nrow=num_states)

# make all diagonal numbers NA
for(i in 1:nrow(best_matrix)){
  best_matrix[i,i] <- 0
}

rownames(best_matrix) <- colnames(best_matrix) <- sort(coding$hab_pref)

# populate matrix with estimated parameters from best diversitree model
best_matrix[best_matrix != 0] <- arrange(diversitree_df, state_2, state_1) %>% pull(transition_rate)

# make diagonal values make things sum to 0
diag(best_matrix) <- -rowSums(best_matrix)

# simmap_best <- make.simmap(tree, hab_pref, nsim = 1000, Q = best_matrix)

# need to split this result up so that files are less than 50MB for GitHub
# number of splits
n_splits <- 10

# find start of each split
splits <- seq(from = 1, to = 1000, by = 1000/10)

# save files out
for(i in 1:length(splits)){
  # save out every 100 sims
  saveRDS(simmap_best[splits[i]:(splits[i]+99)], paste("data/sequencing_rpoB/processed/transition_rates/simmap/simmap_", i, '.rds', sep = ''))
}

# summarise_phytools

## # summarise number of switches between states and time spent in each state
simmap_summary <- describe.simmap(simmap_best, plot=FALSE)

## # remove the tree element - this is the simmap_best
simmap_summary$tree <- NULL

saveRDS(simmap_summary, 'data/sequencing_rpoB/processed/transition_rates/simmap/simmap_summary.rds')

# coerce transitions into dataframe
d_transitions <- as.data.frame(simmap_summary$count, col.names = colnames(simmap_summary$count)) %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(cols = c(everything(), -N, -iter), names_to = 'transition', values_to = 'n_transitions') %>%
  separate(transition, c('state_1', 'state_2'), sep = ',', remove = FALSE) %>%
  mutate(label = gsub(',', ' -> ', transition))

# plot these
group_by(d_transitions, transition) %>%
  filter(max(n_transitions) > 0) %>%
  mutate(average = mean(n_transitions)) %>%
  ungroup() %>%
  mutate(label = fct_reorder(label, desc(average))) %>%
  ggplot() +
  geom_histogram(aes(n_transitions), fill = 'white', col = 'black', binwidth = function(x) (max(x)-min(x))/nclass.FD(x)) +
  facet_wrap(~ label, scales = 'free_x') +
  theme_bw() +
  labs(title = 'Distribution of number of transitions between states',
       subtitle = 'Facets are ordered by common transitions',
       x = 'Number of transitions',
       y = 'Count')


d_transitions_summary <- group_by(d_transitions, iter) %>%
  mutate(prop = n_transitions/sum(n_transitions)) %>%
  group_by(state_1, state_2) %>%
  summarise(ave_num = mean(n_transitions),
            ave_prop = mean(prop),
            lower_ci = quantile(prop, 0.025),
            upper_ci = quantile(prop, 0.975),
            .groups = 'drop') %>%
  arrange(desc(ave_prop)) %>%
  filter(ave_num > 0)

table <- select(d_transitions_summary, state_1, state_2, ave_num, ave_prop) %>%
  mutate(across(starts_with('state'), function(x) gsub('mud_and_shore', 'marine mud', x)),
         across(starts_with('state'), function(x) gsub(':', ' + ', x)),
         ave_prop = round(ave_prop, 2),
         ave_num = round(ave_num, 0))

table_flex <- flextable(table) %>%
  set_header_labels(state_1 = 'from',
                    state_2 = 'to',
                    ave_prop = 'proportion of all transitions',
                    ave_num = 'number of transitions') %>%
  align(align = 'center', part = 'header') %>%
  bold(part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit()

# save out
save_as_image(table_flex, here('plots/sequencing_rpoB/analyses/discrete_character_evolution/proportion_of_transitions.png'), zoom = 3, webshot = 'webshot2')

table_flex

# coerce time into dataframe
d_time <- as.data.frame(simmap_summary$times, col.names = colnames(simmap_summary$times)) %>%
  mutate(n = 1:n()) %>%
  select(-total) %>%
  pivot_longer(cols = c(everything(), -n), names_to = 'state', values_to = 'time') %>%
  group_by(n) %>%
  mutate(prop = time/sum(time)) %>%
  ungroup()

# calculate mean time spent in each state
d_timespent <- group_by(d_time, state) %>%
  summarise(mean = mean(prop), .groups = 'drop')


# next with source sink dynamics using counts
d_source_sink_count <- select(d_transitions_summary, away = state_1, into = state_2, ave_num) %>%
  pivot_longer(cols = c(away, into), names_to = 'direction', values_to = 'habitat_preference') %>%
  group_by(habitat_preference, direction) %>%
  summarise(total_count = sum(ave_num), .groups = 'drop') %>%
  pivot_wider(names_from = direction, values_from = total_count) %>%
  mutate(source_sink1 = away / into,
         source_sink2 = into - away)

table_count <- select(d_source_sink_count, habitat_preference, away, into, source_sink1) %>%
  mutate(across(away:source_sink1, ~round(.x, 2)),
         habitat_preference = gsub(':', ' + ', habitat_preference),
         habitat_preference = gsub('_', ' ', habitat_preference)) %>%
  arrange(desc(source_sink1)) %>%
  flextable(.) %>%
  set_header_labels(habitat_preference = 'habitat preference',
                    source_sink1 = 'source sink ratio') %>%
  align(align = 'center', part = 'all') %>%
  align(align = 'left', part = 'body', j = 1) %>%
  bold(part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit()

save_as_image(table_count, here('plots/sequencing_rpoB/analyses/discrete_character_evolution/source_sink_count.png'), zoom = 3, webshot = 'webshot2')


# can we look at expected numbers of transitions if they were all equally likely? ####

# work out proportion of tips are each habitat preference
d_hab_pref <- group_by(d_meta, habitat_preference) %>%
  tally() %>%
  ungroup() %>%
  mutate(prop = n/sum(n))

# work out expectation
d_expectation <- diversitree_df %>%
  #filter(free_param == 'yes') %>%
  select(state_1, state_2) %>%
  left_join(., select(d_hab_pref, state_1 = habitat_preference, state_1_prop = prop)) %>%
  left_join(., select(d_hab_pref, state_2 = habitat_preference, state_2_prop = prop)) %>%
  mutate(expected_prop = state_1_prop * state_2_prop,
         normalised_expectation = expected_prop/sum(expected_prop))

# observed numbers of transitions
d_expectation <- left_join(d_expectation, select(d_transitions_summary, state_1, state_2, ave_num)) %>%
  mutate(ave_num = replace_na(ave_num, 0),
         tot_transitions = sum(ave_num),
         expected_num = tot_transitions*normalised_expectation,
         ratio = ave_num/expected_num) %>%
  select(state_1, state_2, ave_num, expected_num, ratio)

# make this into a table
table_expect <- mutate(d_expectation,
                       across(ave_num:ratio, ~round(.x, 2)),
                       state_1 = gsub(':', ' + ', state_1),
                       state_1 = gsub('_', ' ', state_1),
                       state_2 = gsub(':', ' + ', state_2),
                       state_2 = gsub('_', ' ', state_2)) %>%
  arrange(desc(ratio)) %>%
  flextable(.) %>%
  set_header_labels(state_1 = 'from',
                    state_2 = 'to',
                    ave_num = 'observed number of transitions',
                    expected_num = 'expected number of transitions',
                    ratio = 'ratio'
  ) %>%
  align(align = 'center', part = 'all') %>%
  align(align = 'left', part = 'body', j = 1) %>%
  bold(part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit()

# save out
save_as_image(table_expect, here('plots/sequencing_rpoB/analyses/discrete_character_evolution/expected_transitions.png'), zoom = 3, webshot = 'webshot2')

#--------------------------------------------------------------#
# are extant generalists generally younger than specialists ####
#--------------------------------------------------------------#

# grab dataset for tree
d_tree <- as_tibble(tree) %>%
  filter(!is.na(label)) %>%
  janitor::clean_names() %>%
  rename(tip_label = label) %>%
  left_join(., d_meta)

ggplot(d_tree, aes(habitat_preference, branch_length)) +
  geom_pretty_boxplot(aes(col = habitat_preference, fill = habitat_preference), show.legend = FALSE) +
  geom_point(shape = 21, fill = 'white', col = 'black', position = position_jitter(width = 0.15), alpha = 0.8) +
  theme_bw(base_size = 16) +
  scale_x_discrete(labels = scales::label_wrap(15)) +
  scale_color_manual(values = cols_hab) +
  scale_fill_manual(values = cols_hab)

# set up correlation matrix for the tree
cor_lambda <- corPagel(value = 1, phy = tree, form = ~tip_label)

d_tree <- tibble(tip_label = tree$tip.label) %>%
  left_join(., d_tree)

# fit phylogenetic generalised linear model
d_tree <- data.frame(d_tree) %>%
  select(tip_label, branch_length, habitat_preference)
row.names(d_tree) <- d_tree$tip_label

d_compare <- caper::comparative.data(phy=tree, data=d_tree, names.col = "tip_label", vcv = TRUE, warn.dropped=TRUE)

mod <- caper::pgls(branch_length ~ habitat_preference, d_compare, lambda = "ML") 

anova(mod)