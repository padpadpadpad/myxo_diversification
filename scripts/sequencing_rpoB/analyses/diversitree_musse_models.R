#---------------------#
# fit musse models ####
#---------------------#

# script to run models of discrete character evolution ####

# load packages ####
librarian::shelf(here, tidyverse, phytools, ggpp, castor, diversitree, hisse)

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
tree <- ape::read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_treepl_cv_node_labels.tre')

# alter tip labels to remove family as they will not link to the distance matrix
# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

tree$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

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

# read in best markov model
best_model <- readRDS('data/sequencing_rpoB/processed/transition_rates/mod_custom_3.rds')

# set up sampling fractions, set them all to 1
sampling_frac <- setNames(rep(0.5, times = 5), sort(unique(hab_pref_num)))

# set up likelihood model for diversitree
lik_musse <- make.musse(tree, hab_pref_num, k = max(hab_pref_num), sampling.f = sampling_frac)

# set constraints for transitions that do not occur
# these are taken from the 0s in best_model
lik_musse <- constrain(lik_musse, 
                       q14~0, q24~0, q25~0, q35~0, q41~0, q42~0, q53~0,
                       q45~0,
                       q54~0)

# we can estimate starting values using starting.point.musse()
start_vals <- starting.point.musse(tree, k = max(hab_pref_num))

# replace the start values with those in the best Mk model
for(i in 1:length(best_model$par)){
  start_vals[names(start_vals) == names(best_model$par)[i]] <- unname(best_model$par[i])
}

# fit musse model
fit_musse <- find.mle(lik_musse, x.init = start_vals[argnames(lik_musse)], method = 'subplex', control = list(maxit = 50000))

# set up NULL models

# remove state dependent speciation and extinction parameters
lik_null_no_sse <- constrain(lik_musse, lambda2 ~ lambda1, lambda3 ~ lambda1, lambda4 ~ lambda1, lambda5 ~ lambda1,
                             mu2 ~ mu1, mu3 ~ mu1, mu4 ~ mu1, mu5~mu1)

# remove only speciation parameters
lik_null_no_ss <- constrain(lik_musse, lambda2 ~ lambda1, lambda3 ~ lambda1, lambda4 ~ lambda1, lambda5 ~ lambda1)

# remove only extinction parameters
lik_null_no_se <- constrain(lik_musse, mu2 ~ mu1, mu3 ~ mu1, mu4 ~ mu1, mu5 ~ mu1)

# fit model
fit_musse_no_sse <- find.mle(lik_null_no_sse, x.init = start_vals[argnames(lik_null_no_sse)], method = 'subplex', control = list(maxit = 100000))
fit_musse_no_ss <- find.mle(lik_null_no_ss, x.init = start_vals[argnames(lik_null_no_ss)],method = 'subplex', control = list(maxit = 50000))
fit_musse_no_se <- find.mle(lik_null_no_se, x.init = start_vals[argnames(lik_null_no_se)], method = 'subplex', control = list(maxit = 100000))

# do model selection

# make AIC table
d_musse_aic <- data.frame(model = c('full_sse', 'no_se', 'no_ss', 'no_sse'), aic = c(AIC(fit_musse), AIC(fit_musse_no_se), AIC(fit_musse_no_ss), AIC(fit_musse_no_sse))) %>%
  dplyr::arrange(., aic) %>%
  mutate(weights = round(MuMIn::Weights(aic), 3))

d_musse_aic

# save out MuSSE models and mcmc
saveRDS(fit_musse, 'data/sequencing_rpoB/processed/transition_rates/asv_musse.rds')
saveRDS(fit_musse_no_se, 'data/sequencing_rpoB/processed/transition_rates/asv_musse_no_se.rds')
saveRDS(fit_musse_no_sse, 'data/sequencing_rpoB/processed/transition_rates/asv_musse_no_sse.rds')
saveRDS(fit_musse_no_ss, 'data/sequencing_rpoB/processed/transition_rates/asv_musse_no_ss.rds')

# read in MuSSE models
fit_musse <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_musse.rds')
fit_musse_no_se <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_musse_no_se.rds')
fit_musse_no_sse <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_musse_no_sse.rds')
fit_musse_no_ss <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_musse_no_ss.rds')

# make AIC table
d_musse_aic <- data.frame(model = c('full_sse', 'no_se', 'no_ss', 'no_sse'), aic = c(AIC(fit_musse), AIC(fit_musse_no_se), AIC(fit_musse_no_ss), AIC(fit_musse_no_sse))) %>%
  dplyr::arrange(., aic) %>%
  mutate(weights = round(MuMIn::Weights(aic), 3))

d_musse_aic
fit_musse_no_se$par
