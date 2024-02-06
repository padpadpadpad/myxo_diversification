# testing for the best transition rate model for the 9 random trees picked from the ASV bootstraps

# what this script does ####

# 1. run the Mk models for the 9 random trees picked from the ASV bootstraps
# 2. create the best model plots as for the 95% and 97.7% data

# load packages ####
librarian::shelf(here, tidyverse, ggtree, ggnewscale, RColorBrewer, patchwork, phytools, diversitree, tidygraph, igraph, ggraph, GGally, ggrepel, flextable, padpadpadpad/MicrobioUoE, cli)

# custom_function ####

source('scripts/sequencing_rpoB/analyses/diversitree_helper_functions.R')

# read in data used throughout

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

# set colours
cols_hab <- readRDS(here('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_colours.rds'))

#----------------#
# bootstrap 1 ####
#----------------#

# read in bootstrapped tree from treepl
tree <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_1_treepl_cv.tre')

ape::ltt.plot(tree, log = 'y')
  
# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

tree$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

# setup for analyses

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

# run diversitree analysis of discrete character evolution

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

# run first four models
mod_er <- find.mle(lik_er, inits_er, method = 'subplex', control = list(maxit = 50000))
mod_ard <- find.mle(lik_ard, inits_ard, method = 'subplex', control = list(maxit = 50000))
mod_trans <- find.mle(lik_transient, inits_trans, method = 'subplex', control = list(maxit = 50000))
mod_sym <- find.mle(lik_sym, inits_sym, method = 'subplex', control = list(maxit = 50000))

# compare models
AIC(mod_sym, mod_ard, mod_er, mod_trans) %>% arrange(AIC)

# sort parameter estimates
mod_ard$par %>% sort()

# get parameters close to 0.
mod_ard$par[mod_ard$par < 1e-03] %>% names(.) %>% sort()

# make custom matrix model
lik_custom1 <- constrain(lik_ard, 
                         q14~0, q24~0, q25~0, q35~0, q41~0, q42~0, q52~0, q54~0)

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
                         q14~0, q24~0, q25~0, q35~0, q41~0, q42~0, q52~0, q54~0,
                         q45~0)

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
                         q14~0, q24~0, q25~0, q35~0, q41~0, q42~0, q52~0, q54~0,
                         q45~0,
                         q43~0)

# make start parameters
inits_custom3 <- rep(1, length(argnames(lik_custom3)))

# run model
mod_custom3 <- find.mle(lik_custom3, inits_custom3, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3) %>% arrange(AIC)

# lets look at model weights
d_aic <- AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3) %>%
  data.frame() %>%
  mutate(log_lik = c(mod_sym$lnLik, mod_ard$lnLik, mod_er$lnLik, mod_trans$lnLik, mod_custom1$lnLik, mod_custom2$lnLik, mod_custom3$lnLik),
         ntips = 1023) %>%
  janitor::clean_names() %>%
  rownames_to_column(var = 'model') %>%
  mutate(aicc = -2*log_lik + 2*df*(ntips/(ntips - df - 1)),
         aic_weight = round(MuMIn::Weights(aic), 2),
         aicc_weight = round(MuMIn::Weights(aicc), 2)) %>%
  arrange(aic)

d_aic

# make this into a table
d_table <- select(d_aic, model, df, log_lik, aic, aic_weight) %>%
  mutate(model = case_when(model == 'mod_er' ~ 'ER',
                           model == 'mod_ard' ~ 'ARD',
                           model == 'mod_sym' ~ 'SYM',
                           model == 'mod_trans' ~ 'Stepwise',
                           model == 'mod_custom1' ~ 'Simplified ARD 1',
                           model == 'mod_custom2' ~ 'Simplified ARD 2',
                           model == 'mod_custom3' ~ 'Simplified ARD 3'))

# make table
table_aic <- d_table %>%
  mutate(across(where(is.numeric), \(x) round(x,2))) %>%
  flextable() %>%
  align(align = 'center', part = 'all') %>%
  set_header_labels(model = "Model",
                    df = 'd.f.',
                    log_lik = 'Log Likelihood',
                    aic = 'AIC',
                    aic_weight = "AIC weight") %>%
  italic(j = 2, part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() 

# make plot of best model
d_habpref_summary <- group_by(d_meta, habitat_preference) %>%
  tally() %>%
  ungroup() %>%
  mutate(prop = n/sum(n))

# make very simple plot to grab legend from - see diversitree_helper_functions.R
p_legend <- make_legend(d_habpref_summary, cols_hab)

# make matrix plot - see diversitree_helper_functions.R
d_diversitree <- get_diversitree_df(mod_custom2, coding$hab_pref_num, coding$hab_pref)

p <- make_network_diversitree(d_diversitree, d_habpref_summary, cols_hab) +
  labs(title = 'Best model from ASV data: random bootstrap 1')

# make source sink table - see diversitree_helper_functions.R
table_rate <- make_source_sink_table(d_diversitree)

# make ltt plot
p_ltt <- plot_ltt(tree, log = 'Y')

# make single plot summarising these results

# save out best model
saveRDS(mod_custom2, 'data/sequencing_rpoB/processed/transition_rates/mod_custom_2_asvboot1.rds')

# make a custom layout for the plot
layout <- c(
  'AAAAB
   AAAAE
   CCC#D'
)

p + 
  p_legend +
  gen_grob(table_aic) + 
  gen_grob(table_rate) +
  p_ltt +
  plot_layout(design = layout, heights = c(0.4, 0.4, 0.2),
              widths = c(0.2, 0.2, 0.2, 0.05, 0.3))

ggsave('plots/manuscript_plots/transitions_asv_boot1.png', height = 6.5, width = 8.5)

#----------------#
# bootstrap 3 ####
#----------------#

# read in bootstrapped tree from treepl
tree <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_3_treepl_cv.tre')

ape::ltt.plot(tree, log = 'y')

# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

tree$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

# setup for analyses
d_meta <- left_join(select(d_habpref, otu, habitat_preference = habitat_preference3, num_present), select(d_taxa, otu:family))

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

# run diversitree analysis of discrete character evolution

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

# run first four models
mod_er <- find.mle(lik_er, inits_er, method = 'subplex', control = list(maxit = 50000))
mod_ard <- find.mle(lik_ard, inits_ard, method = 'subplex', control = list(maxit = 50000))
mod_trans <- find.mle(lik_transient, inits_trans, method = 'subplex', control = list(maxit = 50000))
mod_sym <- find.mle(lik_sym, inits_sym, method = 'subplex', control = list(maxit = 50000))

# compare models
AIC(mod_sym, mod_ard, mod_er, mod_trans) %>% arrange(AIC)

# sort parameter estimates
mod_ard$par %>% sort()

# get parameters close to 0.
mod_ard$par[mod_ard$par < 1e-03] %>% names(.) %>% sort()

# make custom matrix model
lik_custom1 <- constrain(lik_ard, 
                         q14~0, q24~0, q25~0, q31~0, q35~0, q42~0, q45~0, q53~0)

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
                         q14~0, q24~0, q25~0, q31~0, q35~0, q42~0, q45~0, q53~0,
                         q41~0)

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
                         q14~0, q24~0, q25~0, q31~0, q35~0, q42~0, q45~0, q53~0,
                         q41~0,
                         q54~0)

# make start parameters
inits_custom3 <- rep(1, length(argnames(lik_custom3)))

# run model
mod_custom3 <- find.mle(lik_custom3, inits_custom3, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3) %>% arrange(AIC)

# lets look at model weights
d_aic <- AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3) %>%
  data.frame() %>%
  mutate(log_lik = c(mod_sym$lnLik, mod_ard$lnLik, mod_er$lnLik, mod_trans$lnLik, mod_custom1$lnLik, mod_custom2$lnLik, mod_custom3$lnLik),
         ntips = 1023) %>%
  janitor::clean_names() %>%
  rownames_to_column(var = 'model') %>%
  mutate(aicc = -2*log_lik + 2*df*(ntips/(ntips - df - 1)),
         aic_weight = round(MuMIn::Weights(aic), 2),
         aicc_weight = round(MuMIn::Weights(aicc), 2)) %>%
  arrange(aic)

d_aic

# make this into a table
d_table <- select(d_aic, model, df, log_lik, aic, aic_weight) %>%
  mutate(model = case_when(model == 'mod_er' ~ 'ER',
                           model == 'mod_ard' ~ 'ARD',
                           model == 'mod_sym' ~ 'SYM',
                           model == 'mod_trans' ~ 'Stepwise',
                           model == 'mod_custom1' ~ 'Simplified ARD 1',
                           model == 'mod_custom2' ~ 'Simplified ARD 2',
                           model == 'mod_custom3' ~ 'Simplified ARD 3'))

# make table
table_aic <- d_table %>%
  mutate(across(where(is.numeric), \(x) round(x,2))) %>%
  flextable() %>%
  align(align = 'center', part = 'all') %>%
  set_header_labels(model = "Model",
                    df = 'd.f.',
                    log_lik = 'Log Likelihood',
                    aic = 'AIC',
                    aic_weight = "AIC weight") %>%
  italic(j = 2, part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() 

# make plot of best model
d_habpref_summary <- group_by(d_meta, habitat_preference) %>%
  tally() %>%
  ungroup() %>%
  mutate(prop = n/sum(n))

# make very simple plot to grab legend from - see diversitree_helper_functions.R
p_legend <- make_legend(d_habpref_summary, cols_hab)

# make matrix plot - see diversitree_helper_functions.R
d_diversitree <- get_diversitree_df(mod_custom2, coding$hab_pref_num, coding$hab_pref)

p <- make_network_diversitree(d_diversitree, d_habpref_summary, cols_hab) +
  labs(title = 'Best model from ASV data: random bootstrap 3')

# make source sink table - see diversitree_helper_functions.R
table_rate <- make_source_sink_table(d_diversitree)

# make ltt plot
p_ltt <- plot_ltt(tree, log = 'Y')

# make single plot summarising these results

# save out best model
saveRDS(mod_custom2, 'data/sequencing_rpoB/processed/transition_rates/mod_custom_2_asvboot3.rds')

# make a custom layout for the plot
layout <- c(
  'AAAAB
   AAAAE
   CCC#D'
)

p + 
  p_legend +
  gen_grob(table_aic) + 
  gen_grob(table_rate) +
  p_ltt +
  plot_layout(design = layout, heights = c(0.4, 0.4, 0.2),
              widths = c(0.2, 0.2, 0.2, 0.05, 0.3))

ggsave('plots/manuscript_plots/transitions_asv_boot3.png', height = 6.5, width = 8.5)

#----------------#
# bootstrap 4 ####
#----------------#

# read in bootstrapped tree from treepl
tree <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_4_treepl_cv.tre')

ape::ltt.plot(tree, log = 'y')

# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

tree$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

# setup for analyses
d_meta <- left_join(select(d_habpref, otu, habitat_preference = habitat_preference3, num_present), select(d_taxa, otu:family))

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

# run diversitree analysis of discrete character evolution

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

# run first four models
mod_er <- find.mle(lik_er, inits_er, method = 'subplex', control = list(maxit = 50000))
mod_ard <- find.mle(lik_ard, inits_ard, method = 'subplex', control = list(maxit = 50000))
mod_trans <- find.mle(lik_transient, inits_trans, method = 'subplex', control = list(maxit = 50000))
mod_sym <- find.mle(lik_sym, inits_sym, method = 'subplex', control = list(maxit = 50000))

# compare models
AIC(mod_sym, mod_ard, mod_er, mod_trans) %>% arrange(AIC)

# sort parameter estimates
mod_ard$par %>% sort()

# get parameters close to 0.
mod_ard$par[mod_ard$par < 1e-03] %>% names(.) %>% sort()

# make custom matrix model
lik_custom1 <- constrain(lik_ard, 
                         q14~0, q24~0, q25~0, q35~0, q41~0, q42~0, q53~0)

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
                         q14~0, q24~0, q25~0, q35~0, q41~0, q42~0, q53~0,
                         q45~0)

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
                         q14~0, q24~0, q25~0, q35~0, q41~0, q42~0, q53~0,
                         q45~0,
                         q54~0)

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
                         q14~0, q24~0, q25~0, q35~0, q41~0, q42~0, q53~0,
                         q45~0,
                         q54~0,
                         q43~0)

# make start parameters
inits_custom4 <- rep(1, length(argnames(lik_custom4)))

# run model
mod_custom4 <- find.mle(lik_custom4, inits_custom4, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3, mod_custom4) %>% arrange(AIC)

# lets look at model weights
d_aic <- AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3, mod_custom4) %>%
  data.frame() %>%
  mutate(log_lik = c(mod_sym$lnLik, mod_ard$lnLik, mod_er$lnLik, mod_trans$lnLik, mod_custom1$lnLik, mod_custom2$lnLik, mod_custom3$lnLik, mod_custom4$lnLik),
         ntips = 1023) %>%
  janitor::clean_names() %>%
  rownames_to_column(var = 'model') %>%
  mutate(aicc = -2*log_lik + 2*df*(ntips/(ntips - df - 1)),
         aic_weight = round(MuMIn::Weights(aic), 2),
         aicc_weight = round(MuMIn::Weights(aicc), 2)) %>%
  arrange(aic)

d_aic

# make this into a table
d_table <- select(d_aic, model, df, log_lik, aic, aic_weight) %>%
  mutate(model = case_when(model == 'mod_er' ~ 'ER',
                           model == 'mod_ard' ~ 'ARD',
                           model == 'mod_sym' ~ 'SYM',
                           model == 'mod_trans' ~ 'Stepwise',
                           model == 'mod_custom1' ~ 'Simplified ARD 1',
                           model == 'mod_custom2' ~ 'Simplified ARD 2',
                           model == 'mod_custom3' ~ 'Simplified ARD 3',
                           model == 'mod_custom4' ~ 'Simplified ARD 4'))

# make table
table_aic <- d_table %>%
  mutate(across(where(is.numeric), \(x) round(x,2))) %>%
  flextable() %>%
  align(align = 'center', part = 'all') %>%
  set_header_labels(model = "Model",
                    df = 'd.f.',
                    log_lik = 'Log Likelihood',
                    aic = 'AIC',
                    aic_weight = "AIC weight") %>%
  italic(j = 2, part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() 

# make plot of best model
d_habpref_summary <- group_by(d_meta, habitat_preference) %>%
  tally() %>%
  ungroup() %>%
  mutate(prop = n/sum(n))

# make very simple plot to grab legend from - see diversitree_helper_functions.R
p_legend <- make_legend(d_habpref_summary, cols_hab)

# make matrix plot - see diversitree_helper_functions.R
d_diversitree <- get_diversitree_df(mod_custom3, coding$hab_pref_num, coding$hab_pref)

p <- make_network_diversitree(d_diversitree, d_habpref_summary, cols_hab) +
  labs(title = 'Best model from ASV data: random bootstrap 4')

# make source sink table - see diversitree_helper_functions.R
table_rate <- make_source_sink_table(d_diversitree)

# make ltt plot
p_ltt <- plot_ltt(tree, log = 'Y')

# make single plot summarising these results

# save out best model
saveRDS(mod_custom3, 'data/sequencing_rpoB/processed/transition_rates/mod_custom_3_asvboot4.rds')

# make a custom layout for the plot
layout <- c(
  'AAAAB
   AAAAE
   CCC#D'
)

p + 
  p_legend +
  gen_grob(table_aic) + 
  gen_grob(table_rate) +
  p_ltt +
  plot_layout(design = layout, heights = c(0.4, 0.4, 0.2),
              widths = c(0.2, 0.2, 0.2, 0.05, 0.3))

ggsave('plots/manuscript_plots/transitions_asv_boot4.png', height = 6.5, width = 8.5)

#----------------#
# bootstrap 5 ####
#----------------#

# read in bootstrapped tree from treepl
tree <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_5_treepl_cv.tre')

ape::ltt.plot(tree, log = 'y')

# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

tree$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

# setup for analyses
d_meta <- left_join(select(d_habpref, otu, habitat_preference = habitat_preference3, num_present), select(d_taxa, otu:family))

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

# run diversitree analysis of discrete character evolution

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

# run first four models
mod_er <- find.mle(lik_er, inits_er, method = 'subplex', control = list(maxit = 50000))
mod_ard <- find.mle(lik_ard, inits_ard, method = 'subplex', control = list(maxit = 50000))
mod_trans <- find.mle(lik_transient, inits_trans, method = 'subplex', control = list(maxit = 50000))
mod_sym <- find.mle(lik_sym, inits_sym, method = 'subplex', control = list(maxit = 50000))

# compare models
AIC(mod_sym, mod_ard, mod_er, mod_trans) %>% arrange(AIC)

# sort parameter estimates
mod_ard$par %>% sort()

# get parameters close to 0.
mod_ard$par[mod_ard$par < 1e-03] %>% names(.) %>% sort()

# make custom matrix model
lik_custom1 <- constrain(lik_ard, 
                         q14~0, q24~0, q25~0, q31~0, q35~0, q41~0, q42~0, q52~0)

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
                         q14~0, q24~0, q25~0, q31~0, q35~0, q41~0, q42~0, q52~0,
                         q45~0)

# make start parameters
inits_custom2 <- rep(1, length(argnames(lik_custom2)))

# run model
mod_custom2 <- find.mle(lik_custom2, inits_custom2, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2) %>% arrange(AIC)

# lets look at model weights
d_aic <- AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2) %>%
  data.frame() %>%
  mutate(log_lik = c(mod_sym$lnLik, mod_ard$lnLik, mod_er$lnLik, mod_trans$lnLik, mod_custom1$lnLik, mod_custom2$lnLik),
         ntips = 1023) %>%
  janitor::clean_names() %>%
  rownames_to_column(var = 'model') %>%
  mutate(aicc = -2*log_lik + 2*df*(ntips/(ntips - df - 1)),
         aic_weight = round(MuMIn::Weights(aic), 2),
         aicc_weight = round(MuMIn::Weights(aicc), 2)) %>%
  arrange(aic)

d_aic

# make this into a table
d_table <- select(d_aic, model, df, log_lik, aic, aic_weight) %>%
  mutate(model = case_when(model == 'mod_er' ~ 'ER',
                           model == 'mod_ard' ~ 'ARD',
                           model == 'mod_sym' ~ 'SYM',
                           model == 'mod_trans' ~ 'Stepwise',
                           model == 'mod_custom1' ~ 'Simplified ARD 1',
                           model == 'mod_custom2' ~ 'Simplified ARD 2'))

# make table
table_aic <- d_table %>%
  mutate(across(where(is.numeric), \(x) round(x,2))) %>%
  flextable() %>%
  align(align = 'center', part = 'all') %>%
  set_header_labels(model = "Model",
                    df = 'd.f.',
                    log_lik = 'Log Likelihood',
                    aic = 'AIC',
                    aic_weight = "AIC weight") %>%
  italic(j = 2, part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() 

# make plot of best model
d_habpref_summary <- group_by(d_meta, habitat_preference) %>%
  tally() %>%
  ungroup() %>%
  mutate(prop = n/sum(n))

# make very simple plot to grab legend from - see diversitree_helper_functions.R
p_legend <- make_legend(d_habpref_summary, cols_hab)

# make matrix plot - see diversitree_helper_functions.R
d_diversitree <- get_diversitree_df(mod_custom1, coding$hab_pref_num, coding$hab_pref)

p <- make_network_diversitree(d_diversitree, d_habpref_summary, cols_hab) +
  labs(title = 'Best model from ASV data: random bootstrap 5')

# make source sink table - see diversitree_helper_functions.R
table_rate <- make_source_sink_table(d_diversitree)

# make ltt plot
p_ltt <- plot_ltt(tree, log = 'Y')

# make single plot summarising these results

# save out best model
saveRDS(mod_custom1, 'data/sequencing_rpoB/processed/transition_rates/mod_custom_1_asvboot5.rds')

# make a custom layout for the plot
layout <- c(
  'AAAAB
   AAAAE
   CCC#D'
)

p + 
  p_legend +
  gen_grob(table_aic) + 
  gen_grob(table_rate) +
  p_ltt +
  plot_layout(design = layout, heights = c(0.4, 0.4, 0.2),
              widths = c(0.2, 0.2, 0.2, 0.05, 0.3))

ggsave('plots/manuscript_plots/transitions_asv_boot5.png', height = 6.5, width = 8.5)

#----------------#
# bootstrap 6 ####
#----------------#

# read in bootstrapped tree from treepl
tree <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_6_treepl_cv.tre')

ape::ltt.plot(tree, log = 'y')

# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

tree$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

# setup for analyses
d_meta <- left_join(select(d_habpref, otu, habitat_preference = habitat_preference3, num_present), select(d_taxa, otu:family))

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

# run diversitree analysis of discrete character evolution

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

# run first four models
mod_er <- find.mle(lik_er, inits_er, method = 'subplex', control = list(maxit = 50000))
mod_ard <- find.mle(lik_ard, inits_ard, method = 'subplex', control = list(maxit = 50000))
mod_trans <- find.mle(lik_transient, inits_trans, method = 'subplex', control = list(maxit = 50000))
mod_sym <- find.mle(lik_sym, inits_sym, method = 'subplex', control = list(maxit = 50000))

# compare models
AIC(mod_sym, mod_ard, mod_er, mod_trans) %>% arrange(AIC)

# sort parameter estimates
mod_ard$par %>% sort()

# get parameters close to 0.
mod_ard$par[mod_ard$par < 1e-03] %>% names(.) %>% sort()

# make custom matrix model
lik_custom1 <- constrain(lik_ard, 
                         q12~0, q24~0, q25~0, q35~0, q41~0, q42~0, q45~0, q53~0, q54~0)

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
                         q12~0, q24~0, q25~0, q35~0, q41~0, q42~0, q45~0, q53~0, q54~0,
                         q14~0)

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
                         q12~0, q24~0, q25~0, q35~0, q41~0, q42~0, q45~0, q53~0, q54~0,
                         q14~0,
                         q43~0)

# make start parameters
inits_custom3 <- rep(1, length(argnames(lik_custom3)))

# run model
mod_custom3 <- find.mle(lik_custom3, inits_custom3, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3) %>% arrange(AIC)

# lets look at model weights
d_aic <- AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3) %>%
  data.frame() %>%
  mutate(log_lik = c(mod_sym$lnLik, mod_ard$lnLik, mod_er$lnLik, mod_trans$lnLik, mod_custom1$lnLik, mod_custom2$lnLik, mod_custom3$lnLik),
         ntips = 1023) %>%
  janitor::clean_names() %>%
  rownames_to_column(var = 'model') %>%
  mutate(aicc = -2*log_lik + 2*df*(ntips/(ntips - df - 1)),
         aic_weight = round(MuMIn::Weights(aic), 2),
         aicc_weight = round(MuMIn::Weights(aicc), 2)) %>%
  arrange(aic)

d_aic

# make this into a table
d_table <- select(d_aic, model, df, log_lik, aic, aic_weight) %>%
  mutate(model = case_when(model == 'mod_er' ~ 'ER',
                           model == 'mod_ard' ~ 'ARD',
                           model == 'mod_sym' ~ 'SYM',
                           model == 'mod_trans' ~ 'Stepwise',
                           model == 'mod_custom1' ~ 'Simplified ARD 1',
                           model == 'mod_custom2' ~ 'Simplified ARD 2',
                           model == 'mod_custom3' ~ 'Simplified ARD 3'))

# make table
table_aic <- d_table %>%
  mutate(across(where(is.numeric), \(x) round(x,2))) %>%
  flextable() %>%
  align(align = 'center', part = 'all') %>%
  set_header_labels(model = "Model",
                    df = 'd.f.',
                    log_lik = 'Log Likelihood',
                    aic = 'AIC',
                    aic_weight = "AIC weight") %>%
  italic(j = 2, part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() 

# make plot of best model
d_habpref_summary <- group_by(d_meta, habitat_preference) %>%
  tally() %>%
  ungroup() %>%
  mutate(prop = n/sum(n))

# make very simple plot to grab legend from - see diversitree_helper_functions.R
p_legend <- make_legend(d_habpref_summary, cols_hab)

# make matrix plot - see diversitree_helper_functions.R
d_diversitree <- get_diversitree_df(mod_custom2, coding$hab_pref_num, coding$hab_pref)

p <- make_network_diversitree(d_diversitree, d_habpref_summary, cols_hab) +
  labs(title = 'Best model from ASV data: random bootstrap 6')

# make source sink table - see diversitree_helper_functions.R
table_rate <- make_source_sink_table(d_diversitree)

# make ltt plot
p_ltt <- plot_ltt(tree, log = 'Y')

# make single plot summarising these results

# save out best model
saveRDS(mod_custom2, 'data/sequencing_rpoB/processed/transition_rates/mod_custom_2_asvboot6.rds')

# make a custom layout for the plot
layout <- c(
  'AAAAB
   AAAAE
   CCC#D'
)

p + 
  p_legend +
  gen_grob(table_aic) + 
  gen_grob(table_rate) +
  p_ltt +
  plot_layout(design = layout, heights = c(0.4, 0.4, 0.2),
              widths = c(0.2, 0.2, 0.2, 0.05, 0.3))

ggsave('plots/manuscript_plots/transitions_asv_boot6.png', height = 6.5, width = 8.5)

#----------------#
# bootstrap 7 ####
#----------------#

# read in bootstrapped tree from treepl
tree <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_7_treepl_cv.tre')

ape::ltt.plot(tree, log = 'y')

# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

tree$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

# setup for analyses
d_meta <- left_join(select(d_habpref, otu, habitat_preference = habitat_preference3, num_present), select(d_taxa, otu:family))

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

# run diversitree analysis of discrete character evolution

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

# run first four models
mod_er <- find.mle(lik_er, inits_er, method = 'subplex', control = list(maxit = 50000))
mod_ard <- find.mle(lik_ard, inits_ard, method = 'subplex', control = list(maxit = 50000))
mod_trans <- find.mle(lik_transient, inits_trans, method = 'subplex', control = list(maxit = 50000))
mod_sym <- find.mle(lik_sym, inits_sym, method = 'subplex', control = list(maxit = 50000))

# compare models
AIC(mod_sym, mod_ard, mod_er, mod_trans) %>% arrange(AIC)

# sort parameter estimates
mod_ard$par %>% sort()

# get parameters close to 0.
mod_ard$par[mod_ard$par < 1e-03] %>% names(.) %>% sort()

# make custom matrix model
lik_custom1 <- constrain(lik_ard, 
                        q24~0, q25~0, q35~0, q41~0, q42~0, q52~0, q53~0, q54~0)

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
                         q24~0, q25~0, q35~0, q41~0, q42~0, q52~0, q53~0, q54~0,
                         q45~0)

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
                         q14~0, q24~0, q25~0, q31~0, q35~0, q41~0, q42~0, q52~0,
                         q45~0,
                         q14~0)

# make start parameters
inits_custom3 <- rep(1, length(argnames(lik_custom3)))

# run model
mod_custom3 <- find.mle(lik_custom3, inits_custom3, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3) %>% arrange(AIC)

# lets look at model weights
d_aic <- AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3) %>%
  data.frame() %>%
  mutate(log_lik = c(mod_sym$lnLik, mod_ard$lnLik, mod_er$lnLik, mod_trans$lnLik, mod_custom1$lnLik, mod_custom2$lnLik, mod_custom3$lnLik),
         ntips = 1023) %>%
  janitor::clean_names() %>%
  rownames_to_column(var = 'model') %>%
  mutate(aicc = -2*log_lik + 2*df*(ntips/(ntips - df - 1)),
         aic_weight = round(MuMIn::Weights(aic), 2),
         aicc_weight = round(MuMIn::Weights(aicc), 2)) %>%
  arrange(aic)

d_aic

# make this into a table
d_table <- select(d_aic, model, df, log_lik, aic, aic_weight) %>%
  mutate(model = case_when(model == 'mod_er' ~ 'ER',
                           model == 'mod_ard' ~ 'ARD',
                           model == 'mod_sym' ~ 'SYM',
                           model == 'mod_trans' ~ 'Stepwise',
                           model == 'mod_custom1' ~ 'Simplified ARD 1',
                           model == 'mod_custom2' ~ 'Simplified ARD 2',
                           model == 'mod_custom3' ~ 'Simplified ARD 3'))

# make table
table_aic <- d_table %>%
  mutate(across(where(is.numeric), \(x) round(x,2))) %>%
  flextable() %>%
  align(align = 'center', part = 'all') %>%
  set_header_labels(model = "Model",
                    df = 'd.f.',
                    log_lik = 'Log Likelihood',
                    aic = 'AIC',
                    aic_weight = "AIC weight") %>%
  italic(j = 2, part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() 

# make plot of best model
d_habpref_summary <- group_by(d_meta, habitat_preference) %>%
  tally() %>%
  ungroup() %>%
  mutate(prop = n/sum(n))

# make very simple plot to grab legend from - see diversitree_helper_functions.R
p_legend <- make_legend(d_habpref_summary, cols_hab)

# make matrix plot - see diversitree_helper_functions.R
d_diversitree <- get_diversitree_df(mod_custom2, coding$hab_pref_num, coding$hab_pref)

p <- make_network_diversitree(d_diversitree, d_habpref_summary, cols_hab) +
  labs(title = 'Best model from ASV data: random bootstrap 7')

# make source sink table - see diversitree_helper_functions.R
table_rate <- make_source_sink_table(d_diversitree)

# make ltt plot
p_ltt <- plot_ltt(tree, log = 'Y')

# make single plot summarising these results

# save out best model
saveRDS(mod_custom2, 'data/sequencing_rpoB/processed/transition_rates/mod_custom_2_asvboot7.rds')

# make a custom layout for the plot
layout <- c(
  'AAAAB
   AAAAE
   CCC#D'
)

p + 
  p_legend +
  gen_grob(table_aic) + 
  gen_grob(table_rate) +
  p_ltt +
  plot_layout(design = layout, heights = c(0.4, 0.4, 0.2),
              widths = c(0.2, 0.2, 0.2, 0.05, 0.3))

ggsave('plots/manuscript_plots/transitions_asv_boot7.png', height = 6.5, width = 8.5)

#----------------#
# bootstrap 8 ####
#----------------#

# read in bootstrapped tree from treepl
tree <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_8_treepl_cv.tre')

ape::ltt.plot(tree, log = 'y')

# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

tree$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

# setup for analyses
d_meta <- left_join(select(d_habpref, otu, habitat_preference = habitat_preference3, num_present), select(d_taxa, otu:family))

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

# run diversitree analysis of discrete character evolution

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

# run first four models
mod_er <- find.mle(lik_er, inits_er, method = 'subplex', control = list(maxit = 50000))
mod_ard <- find.mle(lik_ard, inits_ard, method = 'subplex', control = list(maxit = 50000))
mod_trans <- find.mle(lik_transient, inits_trans, method = 'subplex', control = list(maxit = 50000))
mod_sym <- find.mle(lik_sym, inits_sym, method = 'subplex', control = list(maxit = 50000))

# compare models
AIC(mod_sym, mod_ard, mod_er, mod_trans) %>% arrange(AIC)

# sort parameter estimates
mod_ard$par %>% sort()

# get parameters close to 0.
mod_ard$par[mod_ard$par < 1e-03] %>% names(.) %>% sort()

# make custom matrix model
lik_custom1 <- constrain(lik_ard, 
                         q14~0, q24~0, q25~0, q31~0, q35~0, q41~0, q42~0, q52~0, q53~0)

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
                         q14~0, q24~0, q25~0, q31~0, q35~0, q41~0, q42~0, q52~0, q53~0,
                         q54~0)

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
                         q14~0, q24~0, q25~0, q31~0, q35~0, q41~0, q42~0, q52~0, q53~0,
                         q54~0,
                         q45~0)

# make start parameters
inits_custom3 <- rep(1, length(argnames(lik_custom3)))

# run model
mod_custom3 <- find.mle(lik_custom3, inits_custom3, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3) %>% arrange(AIC)

# lets look at model weights
d_aic <- AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3) %>%
  data.frame() %>%
  mutate(log_lik = c(mod_sym$lnLik, mod_ard$lnLik, mod_er$lnLik, mod_trans$lnLik, mod_custom1$lnLik, mod_custom2$lnLik, mod_custom3$lnLik),
         ntips = 1023) %>%
  janitor::clean_names() %>%
  rownames_to_column(var = 'model') %>%
  mutate(aicc = -2*log_lik + 2*df*(ntips/(ntips - df - 1)),
         aic_weight = round(MuMIn::Weights(aic), 2),
         aicc_weight = round(MuMIn::Weights(aicc), 2)) %>%
  arrange(aic)

d_aic

# make this into a table
d_table <- select(d_aic, model, df, log_lik, aic, aic_weight) %>%
  mutate(model = case_when(model == 'mod_er' ~ 'ER',
                           model == 'mod_ard' ~ 'ARD',
                           model == 'mod_sym' ~ 'SYM',
                           model == 'mod_trans' ~ 'Stepwise',
                           model == 'mod_custom1' ~ 'Simplified ARD 1',
                           model == 'mod_custom2' ~ 'Simplified ARD 2',
                           model == 'mod_custom3' ~ 'Simplified ARD 3'))

# make table
table_aic <- d_table %>%
  mutate(across(where(is.numeric), \(x) round(x,2))) %>%
  flextable() %>%
  align(align = 'center', part = 'all') %>%
  set_header_labels(model = "Model",
                    df = 'd.f.',
                    log_lik = 'Log Likelihood',
                    aic = 'AIC',
                    aic_weight = "AIC weight") %>%
  italic(j = 2, part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() 

# make plot of best model
d_habpref_summary <- group_by(d_meta, habitat_preference) %>%
  tally() %>%
  ungroup() %>%
  mutate(prop = n/sum(n))

# make very simple plot to grab legend from - see diversitree_helper_functions.R
p_legend <- make_legend(d_habpref_summary, cols_hab)

# make matrix plot - see diversitree_helper_functions.R
d_diversitree <- get_diversitree_df(mod_custom2, coding$hab_pref_num, coding$hab_pref)

p <- make_network_diversitree(d_diversitree, d_habpref_summary, cols_hab) +
  labs(title = 'Best model from ASV data: random bootstrap 8')

# make source sink table - see diversitree_helper_functions.R
table_rate <- make_source_sink_table(d_diversitree)

# make ltt plot
p_ltt <- plot_ltt(tree, log = 'Y')

# make single plot summarising these results

# save out best model
saveRDS(mod_custom2, 'data/sequencing_rpoB/processed/transition_rates/mod_custom_2_asvboot8.rds')

# make a custom layout for the plot
layout <- c(
  'AAAAB
   AAAAE
   CCC#D'
)

p + 
  p_legend +
  gen_grob(table_aic) + 
  gen_grob(table_rate) +
  p_ltt +
  plot_layout(design = layout, heights = c(0.4, 0.4, 0.2),
              widths = c(0.2, 0.2, 0.2, 0.05, 0.3))

ggsave('plots/manuscript_plots/transitions_asv_boot8.png', height = 6.5, width = 8.5)

#----------------#
# bootstrap 9 ####
#----------------#

# read in bootstrapped tree from treepl
tree <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps/raxml_9_treepl_cv.tre')

ape::ltt.plot(tree, log = 'y')

# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

tree$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

# setup for analyses
d_meta <- left_join(select(d_habpref, otu, habitat_preference = habitat_preference3, num_present), select(d_taxa, otu:family))

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

# run diversitree analysis of discrete character evolution

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

# run first four models
mod_er <- find.mle(lik_er, inits_er, method = 'subplex', control = list(maxit = 50000))
mod_ard <- find.mle(lik_ard, inits_ard, method = 'subplex', control = list(maxit = 50000))
mod_trans <- find.mle(lik_transient, inits_trans, method = 'subplex', control = list(maxit = 50000))
mod_sym <- find.mle(lik_sym, inits_sym, method = 'subplex', control = list(maxit = 50000))

# compare models
AIC(mod_sym, mod_ard, mod_er, mod_trans) %>% arrange(AIC)

# sort parameter estimates
mod_ard$par %>% sort()

# get parameters close to 0.
mod_ard$par[mod_ard$par < 1e-03] %>% names(.) %>% sort()

# make custom matrix model
lik_custom1 <- constrain(lik_ard, 
                         q12~0, q13~0, q14~0, q25~0, q31~0, q34~0, q35~0, q41~0, q42~0, q54~0)

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
                         q12~0, q13~0, q14~0, q25~0, q31~0, q34~0, q35~0, q41~0, q42~0, q54~0,
                         q45~0)

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
                         q12~0, q13~0, q14~0, q25~0, q31~0, q34~0, q35~0, q41~0, q42~0, q54~0,
                         q45~0,
                         q43~0)

# make start parameters
inits_custom3 <- rep(1, length(argnames(lik_custom3)))

# run model
mod_custom3 <- find.mle(lik_custom3, inits_custom3, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3) %>% arrange(AIC)

# lets look at model weights
d_aic <- AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3) %>%
  data.frame() %>%
  mutate(log_lik = c(mod_sym$lnLik, mod_ard$lnLik, mod_er$lnLik, mod_trans$lnLik, mod_custom1$lnLik, mod_custom2$lnLik, mod_custom3$lnLik),
         ntips = 1023) %>%
  janitor::clean_names() %>%
  rownames_to_column(var = 'model') %>%
  mutate(aicc = -2*log_lik + 2*df*(ntips/(ntips - df - 1)),
         aic_weight = round(MuMIn::Weights(aic), 2),
         aicc_weight = round(MuMIn::Weights(aicc), 2)) %>%
  arrange(aic)

d_aic

# make this into a table
d_table <- select(d_aic, model, df, log_lik, aic, aic_weight) %>%
  mutate(model = case_when(model == 'mod_er' ~ 'ER',
                           model == 'mod_ard' ~ 'ARD',
                           model == 'mod_sym' ~ 'SYM',
                           model == 'mod_trans' ~ 'Stepwise',
                           model == 'mod_custom1' ~ 'Simplified ARD 1',
                           model == 'mod_custom2' ~ 'Simplified ARD 2',
                           model == 'mod_custom3' ~ 'Simplified ARD 3'))

# make table
table_aic <- d_table %>%
  mutate(across(where(is.numeric), \(x) round(x,2))) %>%
  flextable() %>%
  align(align = 'center', part = 'all') %>%
  set_header_labels(model = "Model",
                    df = 'd.f.',
                    log_lik = 'Log Likelihood',
                    aic = 'AIC',
                    aic_weight = "AIC weight") %>%
  italic(j = 2, part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() 

# make plot of best model
d_habpref_summary <- group_by(d_meta, habitat_preference) %>%
  tally() %>%
  ungroup() %>%
  mutate(prop = n/sum(n))

# make very simple plot to grab legend from - see diversitree_helper_functions.R
p_legend <- make_legend(d_habpref_summary, cols_hab)

# make matrix plot - see diversitree_helper_functions.R
d_diversitree <- get_diversitree_df(mod_custom2, coding$hab_pref_num, coding$hab_pref)

p <- make_network_diversitree(d_diversitree, d_habpref_summary, cols_hab) +
  labs(title = 'Best model from ASV data: random bootstrap 9')

# make source sink table - see diversitree_helper_functions.R
table_rate <- make_source_sink_table(d_diversitree)

# make ltt plot
p_ltt <- plot_ltt(tree, log = 'Y')

# make single plot summarising these results

# save out best model
saveRDS(mod_custom2, 'data/sequencing_rpoB/processed/transition_rates/mod_custom_2_asvboot9.rds')

# make a custom layout for the plot
layout <- c(
  'AAAAB
   AAAAE
   CCC#D'
)

p + 
  p_legend +
  gen_grob(table_aic) + 
  gen_grob(table_rate) +
  p_ltt +
  plot_layout(design = layout, heights = c(0.4, 0.4, 0.2),
              widths = c(0.2, 0.2, 0.2, 0.05, 0.3))

ggsave('plots/manuscript_plots/transitions_asv_boot9.png', height = 6.5, width = 8.5)
