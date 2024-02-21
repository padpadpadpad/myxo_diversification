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
tree <- ape::read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_treepl_cv_node_labels.tre')

# alter tip labels to remove family as they will not link to the distance matrix
# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

tree$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

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
sampling_frac <- setNames(rep(0.5, times = 5), sort(unique(hab_pref_num)))

# read in best Mk model
fit_mk <- readRDS('data/sequencing_rpoB/processed/transition_rates/mod_custom_3.rds')

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
                       q14~0, q24~0, q25~0, q35~0, q41~0, q42~0, q53~0,
                       q45~0,
                       q54~0)

# set all q parameters to be that from the Mk model
lik_musse2 <- constrain(lik_musse,
                        q12~1.2742403,
                        q13~2.1371199,
                        q15~25.3769010,
                        q21~1.7189816,
                        q23~3.9420235,
                        q31~7.4760225,
                        q32~20.7244840,
                        q34~4.5038291,
                        q43~0.2907249,
                        q51~16.4972949,
                        q52~1.0257041)

fit_musse_no_se2 <- find.mle(lik_musse2, x.init = start_vals[argnames(lik_musse2)], method = 'subplex', control = list(maxit = 100000))

# set all q parameters to be the same from the Mk model
lik_musse3 <- constrain(lik_musse,
                        q13~q12,
                        q15~q12,
                        q21~q12,
                        q23~q12,
                        q31~q12,
                        q32~q12,
                        q34~q12,
                        q43~q12,
                        q51~q12,
                        q52~q12)

fit_musse_no_se3 <- find.mle(lik_musse3, x.init = start_vals[argnames(lik_musse3)], method = 'subplex', control = list(maxit = 100000))

data.frame(aic = c(AIC(fit_musse_no_se), AIC(fit_musse_no_se2), AIC(fit_musse_no_se3))) %>%
  mutate(weights = round(MuMIn::Weights(aic), 3))

fit_musse_no_se$par %>% round(2)
fit_musse_no_se2$par %>% round(2)
fit_musse_no_se3$par %>% round(2)

fit_musse_no_se2$lnLik

# compare values 
d_lambda <- data.frame(original = fit_musse_no_se$par[str_detect(names(fit_musse_no_se$par), 'lambda')],
                       new = fit_musse_no_se2$par[str_detect(names(fit_musse_no_se$par), 'lambda')])

d_lambda = mutate(d_lambda, diff = abs(original - new),
                  diff/original * 100)

cor(d_lambda$original, d_lambda$new)

d_q <- data.frame(markov = fit_mk$par[str_detect(names(fit_mk$par), 'q')],
                  musse = fit_musse_no_se$par[str_detect(names(fit_musse_no_se$par), 'q')])
cor(d_q$markov, d_q$musse)

# 