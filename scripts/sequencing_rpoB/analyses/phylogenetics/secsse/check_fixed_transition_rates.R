# for the MuSSE models, compare transition rates to the best Mk model

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
tree <- read.tree(here('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_chronopl10.tre'))

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

# set up sampling fractions, set them all to 1
sampling_frac <- setNames(rep(1, times = 5), sort(unique(hab_pref_num)))

# read in best Mk model
fit_mk <- readRDS('data/sequencing_rpoB/processed/transition_rates/mod_custom5.rds')

fit_musse_no_se <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_musse_no_se.rds')

fit_mk$par
fit_musse_no_se$par
fit_musse_no_se$par[names(fit_musse_no_se$par) %in% names(fit_mk$par)]

d <- tibble(mk = fit_mk$par,
            musse = fit_musse_no_se$par[names(fit_musse_no_se$par) %in% names(fit_mk$par)],
            transition = names(fit_mk$par))

ggplot(d, aes(mk, musse)) +
  geom_point() +
  geom_abline(aes(slope = 1, intercept = 0))

# set up likelihood model for diversitree
lik_musse <- make.musse(tree, hab_pref_num, k = max(hab_pref_num), sampling.f = sampling_frac)

# we can estimate starting values using starting.point.musse()
start_vals <- starting.point.musse(tree, k = max(hab_pref_num))

# set constraints for transitions that do not occur
# these are taken from the 0s in best_model
lik_musse <- constrain(lik_musse,
                        mu2 ~ mu1, mu3 ~ mu1, mu4 ~ mu1, mu5 ~ mu1,
                        q14~0, q24~0, q25~0, q41~0, q53~0,
                        q42~0,
                        q45~0,
                        q54~0,
                        q52~0)

# set all q parameters to be that from the Mk model
lik_musse2 <- constrain(lik_musse,
                        q12~1.4452899,
                        q13~1.4372485,
                        q15~12.4165517,
                        q21~0.3867756,
                        q23~1.6808420,
                        q31~5.0890002,
                        q32~8.2505921,
                        q34~2.4666607,
                        q35~1.7432721,
                        q43~0.2259779,
                        q51~8.9171032)

fit_musse_no_se2 <- find.mle(lik_musse2, x.init = start_vals[argnames(lik_musse2)], method = 'subplex', control = list(maxit = 100000))

# set all q parameters to be that from the Mk model
lik_musse3 <- constrain(lik_musse,
                        q13~q12,
                        q15~q12,
                        q21~q12,
                        q23~q12,
                        q31~q12,
                        q32~q12,
                        q34~q12,
                        q35~q12,
                        q43~q12,
                        q51~q12)

fit_musse_no_se3 <- find.mle(lik_musse3, x.init = start_vals[argnames(lik_musse3)], method = 'subplex', control = list(maxit = 100000))

data.frame(aic = c(AIC(fit_musse_no_se), AIC(fit_musse_no_se2), AIC(fit_musse_no_se3))) %>%
  dplyr::arrange(., aic) %>%
  mutate(weights = round(MuMIn::Weights(aic), 3))

fit_musse_no_se$par %>% round(2)
fit_musse_no_se2$par %>% round(2)
fit_musse_no_se3$par %>% round(2)

fit_musse_no_se2$lnLik
