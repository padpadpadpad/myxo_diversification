
# load packages
library(here)
library(tidyverse)
library(ggtree)
library(ggnewscale)
library(RColorBrewer)
library(patchwork)
library(phytools)
library(MetBrewer)
library(ggpp)
library(castor)
library(diversitree)
library(tidygraph)
library(igraph)
library(ggraph)
library(GGally)
library(ggrepel)
library(flextable)
library(ggridges)
library(hisse)

# set where I am in the project
here::i_am('scripts/sequencing_rpoB/analyses/phylogenetics/estimating_transition_rates_v3.qmd')

# read in habitat preference
d_habpref <- read.csv(here('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_asv.csv'))

# read in phyloseq object and grab tax table
d_taxa <- readRDS(here('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_asv_filt.rds')) %>%
  phyloseq::tax_table() %>%
  data.frame() %>%
  janitor::clean_names() %>%
  rownames_to_column('otu')

# create d_meta
d_meta <- left_join(select(d_habpref, otu, habitat_preference, num_present), select(d_taxa, otu:family))

# read in tree
tree <- read.tree(here('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_chronopl10.tre'))

# read in shift nodes
shiftnodes <- readRDS(here('data/sequencing_rpoB/processed/shiftnodes.rds'))

## ----custom_function------------------------------------------------------------------------------------------------------------

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

## ----tip_state_dataframe--------------------------------------------------------------------------------------------------------
d_meta <- tibble(tip_label = tree$tip.label) %>%
  left_join(., rename(d_meta, tip_label = otu))

# make sure order of habitat preference is the same in the trait vector as the tip labels
sum(d_meta$tip_label == tree$tip.label) == length(tree$tip.label)
# SUCCESS if TRUE

# replace NA of outgroup - make it freshwater:terrestrial - the most common state
# phytools::make.simmap cannot take NAs
hab_pref <- setNames(d_meta$habitat_preference, d_meta$tip_label)
hab_pref[is.na(hab_pref)] <- 'freshwater:terrestrial'

# make habitat preference vector numeric
hab_pref_num <- as.numeric(as.factor(hab_pref)) -1
hab_pref_num <- setNames(hab_pref_num, d_meta$tip_label)
hab_pref_num2 <- hab_pref_num + 1

# coding from numeric to character
coding <- tibble(hab_pref = unname(hab_pref), hab_pref_num = unname(hab_pref_num)) %>%
  distinct() %>%
  mutate(hab_pref2 = gsub(':', '.', hab_pref),
         hab_pref_num2 = hab_pref_num + 1) %>%
  arrange(hab_pref) %>%
  mutate(hab_pref_axis = gsub(':', '/ ', hab_pref),
  hab_pref_axis = gsub('_', ' ', hab_pref_axis),
  # rename the columns for easy naming
  initials = c('F', 'FM', 'FT', 'G', 'M', 'MT', 'T'))

coding


## ----diversitree_load-----------------------------------------------------------------------------------------------------------
# all of these methods needs a likelihood function, we can build a Mkn model
lik_ard <- make.mkn(tree, hab_pref_num2, k = max(hab_pref_num2))

argnames(lik_ard)

# make symmetric rates model
lik_sym <- constrain(lik_ard, 
                     q12~q21, q13~q31, q14~q41, q15~q51, q16~q61, q17~q71,
                     q23~q32, q24~q42, q25~q52, q26~q62, q27~q72,
                     q34~q43, q35~q53, q36~q63, q37~q73,
                     q45~q54, q46~q64, q47~q74,
                     q56~q65, q57~q75,
                     q67~q76
                     )

argnames(lik_sym)

# make equal rates model
lik_er <- constrain.i(lik_ard, rep('q12', length(argnames(lik_ard))), i.free = 1)

lik_er <- constrain(lik_ard,
                    q13~q12, q14~q12, q15~q12, q16~q12, q17~q12,
                    q21~q12, q23~q12, q24~q12, q25~q12, q26~q12, q27~q12,
                    q31~q12, q32~q12, q34~q12, q35~q12, q36~q12, q37~q12,
                    q41~q12, q42~q12, q43~q12, q45~q12, q46~q12, q47~q12,
                    q51~q12, q52~q12, q53~q12, q54~q12, q56~q12, q57~q12,
                    q61~q12, q62~q12, q63~q12, q64~q12, q65~q12, q67~q12,
                    q71~q12, q72~q12, q73~q12, q74~q12, q75~q12, q76~q12)

argnames(lik_er)

# need to pass start values to it - can grab these from the ape::ace, but we will just pass an average rate to the model
inits_ard <- rep(1, length(argnames(lik_ard)))
inits_sym <- rep(1, length(argnames(lik_sym)))
inits_er <- rep(1, length(argnames(lik_er)))

## ----diversitree_compare--------------------------------------------------------------------------------------------------------

AIC(mod_sym, mod_ard, mod_er) %>% arrange(AIC)

anova(mod_ard, mod_sym)
anova(mod_ard, mod_er)



## ----plot_ard-------------------------------------------------------------------------------------------------------------------
# look at paraemters
mod_ard$par

# plot of all parameters
p1 <- tibble(rate = mod_ard$par) %>%
  ggplot(., aes(rate)) +
  geom_histogram(fill = 'white', col = 'black') +
  theme_bw()

# plot of constrained x axis (rates <1)
p2 <- tibble(rate = mod_ard$par) %>%
  filter(rate < 1) %>%
  ggplot(., aes(rate)) +
  geom_histogram(fill = 'white', col = 'black') +
  theme_bw()

p1+p2

# sort parameter estimates
mod_ard$par %>% sort()

# get parameters close to 0.
mod_ard$par[mod_ard$par < 1e-03] %>% names(.)

# make custom matrix model
lik_custom1 <- constrain(lik_ard, 
                     q13~0, q15~0, q16~0, q17~0, q26~0, q45~0, q51~0, q53~0, q54~0, q62~0)

# make start parameters
inits_custom1 <- rep(1, length(argnames(lik_custom1)))

mod_custom1 <- find.mle(lik_custom1, inits_custom1, method = 'subplex', control = list(maxit = 50000))


## ----diversitree_1_compare------------------------------------------------------------------------------------------------------
# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_custom1) %>% arrange(AIC)

# anova
anova(mod_ard, mod_custom1)


## ----plot_custom1---------------------------------------------------------------------------------------------------------------
# filter for only estimated parameters
mod_custom1$par.full[names(mod_custom1$par.full) %in% argnames(lik_custom1)] 

# plot of all parameters
p1 <- tibble(rate = mod_custom1$par.full[names(mod_custom1$par.full) %in% argnames(lik_custom1)]) %>%
  ggplot(., aes(rate)) +
  geom_histogram(fill = 'white', col = 'black') +
  theme_bw()

# plot of constrained x axis (rates <1)
p2 <- tibble(rate = mod_custom1$par.full[names(mod_custom1$par.full) %in% argnames(lik_custom1)]) %>%
  filter(rate < 1) %>%
  ggplot(., aes(rate)) +
  geom_histogram(fill = 'white', col = 'black') +
  theme_bw()

p1+p2

# sort parameter estimates
mod_custom1$par %>% sort()

# get parameters close to 0.
mod_custom1$par[mod_custom1$par < 1e-03] %>% names(.)


## ----diversitree_2_load---------------------------------------------------------------------------------------------------------
# make custom matrix model
lik_custom2 <- constrain(lik_ard, 
                     q13~0, q15~0, q16~0, q17~0, q26~0, q45~0, q51~0, q53~0, q54~0, q62~0,
                     q31~0, q35~0, q36~0, q63~0, q72~0, q73~0, q75~0)

# make start parameters
inits_custom2 <- rep(1, length(argnames(lik_custom2)))

## # run model
mod_custom2 <- find.mle(lik_custom2, inits_custom2, method = 'subplex', control = list(maxit = 50000))


## ----diversitree_2_compare------------------------------------------------------------------------------------------------------
# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_custom1, mod_custom2) %>% arrange(AIC)

# anova
anova(mod_custom1, mod_custom2)

# sort parameter estimates
mod_custom2$par %>% sort()


## ----diversitree_3_load---------------------------------------------------------------------------------------------------------
# make custom matrix model
lik_custom3 <- constrain(lik_ard, 
                     q13~0, q15~0, q16~0, q17~0, q26~0, q45~0, q51~0, q53~0, q54~0, q62~0,
                     q31~0, q35~0, q36~0, q63~0, q72~0, q73~0, q75~0,
                     q57~0)

# make start parameters
inits_custom3 <- rep(1, length(argnames(lik_custom3)))

# run model
mod_custom3 <- find.mle(lik_custom3, inits_custom3, method = 'subplex', control = list(maxit = 50000))

## ----diversitree_3_compare------------------------------------------------------------------------------------------------------
# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_custom1, mod_custom2, mod_custom3) %>% arrange(AIC)

# anova
anova(mod_custom2, mod_custom3)

# sort parameter estimates
mod_custom3$par %>% sort()


## ----diversitree_4_load---------------------------------------------------------------------------------------------------------
# make custom matrix model
lik_custom4 <- constrain(lik_ard, 
                     q13~0, q15~0, q16~0, q17~0, q26~0, q45~0, q51~0, q53~0, q54~0, q62~0,
                     q31~0, q35~0, q36~0, q63~0, q72~0, q73~0, q75~0,
                     q57~0,
                     q43~0)

# make start parameters
inits_custom4 <- rep(1, length(argnames(lik_custom4)))

# run model
mod_custom4 <- find.mle(lik_custom4, inits_custom4, method = 'subplex', control = list(maxit = 50000))

## ----diversitree_4_compare------------------------------------------------------------------------------------------------------
# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_custom1, mod_custom2, mod_custom3, mod_custom4) %>% arrange(AIC)

# anova
anova(mod_custom3, mod_custom4)

# sort parameter estimates
mod_custom4$par %>% sort()


## ----diversitree_5_load---------------------------------------------------------------------------------------------------------
# make custom matrix model
lik_custom5 <- constrain(lik_ard, 
                     q13~0, q15~0, q16~0, q17~0, q26~0, q45~0, q51~0, q53~0, q54~0, q62~0,
                     q31~0, q35~0, q36~0, q63~0, q72~0, q73~0, q75~0,
                     q57~0,
                     q43~0,
                     q56~0)

# make start parameters
inits_custom5 <- rep(1, length(argnames(lik_custom5)))

# run model
mod_custom5 <- find.mle(lik_custom5, inits_custom5, method = 'subplex', control = list(maxit = 50000))

## ----diversitree_5_compare------------------------------------------------------------------------------------------------------
# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_custom1, mod_custom2, mod_custom3, mod_custom4, mod_custom5) %>% arrange(AIC)

# anova
anova(mod_custom4, mod_custom5)

# sort parameter estimates
mod_custom5$par %>% sort()


## ----diversitree_6_load---------------------------------------------------------------------------------------------------------
# read in files
mod_custom6 <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/mod_custom6_mac.rds'))

# make custom matrix model
lik_custom6 <- constrain(lik_ard, 
                     q13~0, q15~0, q16~0, q17~0, q26~0, q45~0, q51~0, q53~0, q54~0, q62~0,
                     q31~0, q35~0, q36~0, q63~0, q72~0, q73~0, q75~0,
                     q57~0,
                     q43~0,
                     q56~0,
                     q27~0)

# make start parameters
inits_custom6 <- rep(1, length(argnames(lik_custom6)))

#-----------------------#
# STOP AT THIS POINT ####
#-----------------------#

# run model - got stuck here
mod_custom6 <- find.mle(lik_custom6, inits_custom6, method = 'subplex', control = list(maxit = 5))

## ----diversitree_6_compare------------------------------------------------------------------------------------------------------
# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_custom1, mod_custom2, mod_custom3, mod_custom4, mod_custom5, mod_custom6) %>% arrange(AIC)

# anova
anova(mod_custom5, mod_custom6)

# sort parameter estimates
mod_custom6$par %>% sort()


## ----diversitree_7_load---------------------------------------------------------------------------------------------------------
# read in files
mod_custom7 <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/mod_custom7_mac.rds'))

# make custom matrix model
lik_custom7 <- constrain(lik_ard, 
                     q16~0, q26~0, q27~0, q35~0, q36~0, q42~0, q46~0, q47~0, q53~0, q54~0, q63~0, q64~0, q72~0, q74~0, q17~0, q75~0, q71~0, q51~0, q14~0, q52~0, q67~0)

# make start parameters
inits_custom7 <- rep(1, length(argnames(lik_custom7)))


## ----diversitree_7_setup--------------------------------------------------------------------------------------------------------
## # make custom matrix model
## lik_custom7 <- constrain(lik_ard,
##                      q16~0, q26~0, q27~0, q35~0, q36~0, q42~0, q46~0, q47~0, q53~0, q54~0, q63~0, q64~0, q72~0, q74~0, q17~0, q75~0, q71~0, q51~0, q14~0, q52~0, q67~0)
## 
## # make start parameters
## inits_custom7 <- rep(1, length(argnames(lik_custom7)))
## 
## # run model
## mod_custom7 <- find.mle(lik_custom7, inits_custom7, method = 'subplex', control = list(maxit = 50000))


## ----diversitree_7_compare------------------------------------------------------------------------------------------------------
# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_custom1, mod_custom2, mod_custom3, mod_custom4, mod_custom5, mod_custom6, mod_custom7) %>% arrange(AIC)

# anova
anova(mod_custom6, mod_custom7)

# sort parameter estimates
mod_custom7$par %>% sort()


## ----diversitree_8_load---------------------------------------------------------------------------------------------------------
# read in files
mod_custom8 <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/mod_custom8_mac.rds'))

# make custom matrix model
lik_custom8 <- constrain(lik_ard, 
                     q16~0, q26~0, q27~0, q35~0, q36~0, q42~0, q46~0, q47~0, q53~0, q54~0, q63~0, q64~0, q72~0, q74~0, q17~0, q75~0, q71~0, q51~0, q14~0, q52~0, q67~0, q34~0)

# make start parameters
inits_custom8 <- rep(1, length(argnames(lik_custom8)))


## ----diversitree_8_setup--------------------------------------------------------------------------------------------------------
## # make custom matrix model
## lik_custom8 <- constrain(lik_ard,
##                      q16~0, q26~0, q27~0, q35~0, q36~0, q42~0, q46~0, q47~0, q53~0, q54~0, q63~0, q64~0, q72~0, q74~0, q17~0, q75~0, q71~0, q51~0, q14~0, q52~0, q67~0, q34~0)
## 
## # make start parameters
## inits_custom8 <- rep(1, length(argnames(lik_custom8)))
## 
## # run model
## mod_custom8 <- find.mle(lik_custom8, inits_custom8, method = 'subplex', control = list(maxit = 50000))


## ----diversitree_8_compare------------------------------------------------------------------------------------------------------
# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_custom1, mod_custom2, mod_custom3, mod_custom4, mod_custom5, mod_custom6, mod_custom7, mod_custom8) %>% arrange(AIC)

# anova
anova(mod_custom7, mod_custom8)

# sort parameter estimates
mod_custom8$par %>% sort()


## ----plot_transition_matrix-----------------------------------------------------------------------------------------------------
diversitree_df <- get_diversitree_df(mod_custom7, coding$hab_pref_num2, coding$hab_pref)

diversitree_df %>%
  left_join(., select(coding, state_1 = hab_pref, state_1_num = hab_pref_num2, state_1_label = hab_pref_axis)) %>%
  left_join(., select(coding, state_2 = hab_pref, state_2_num = hab_pref_num2, state_2_label = hab_pref_axis)) %>%
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
  labs(y = 'state 1',
  x = 'state 2',
  title = paste('all rates different with', length(mod_custom7$par), 'free parameters', sep = ' ')) +
  coord_fixed() +
  scale_color_manual(values = c('red', 'black'))

ggsave(here('plots/sequencing_rpoB/analyses/transition_matrix_diversitree.png'), last_plot(), height = 6, width = 8)



## ----mcmc_load------------------------------------------------------------------------------------------------------------------
fit_mcmc3 <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/best_diversitree_mcmc.rds'))


## ----mcmc_try-------------------------------------------------------------------------------------------------------------------
## # set up initial start values
## inits_mcmc <- mod_custom7$par
## 
## # set up upper and lower limits - limit the values to be < 10 times the max value
## lower_mcmc <- rep(0, length(inits_mcmc))
## upper_mcmc <- rep(max(inits_mcmc)*10, length(inits_mcmc))
## 
## # run first mcmc to tune w
## fit_mcmc <- mcmc(lik_custom7, inits_mcmc, nsteps = 10, w = 0.1, upper = upper_mcmc, lower = lower_mcmc)
## 
## # tune w for each parameter
## w <- diff(sapply(fit_mcmc[2:(ncol(fit_mcmc)-1)], quantile, c(.05, .95)))
## 
## # run second mcmc to tune w
## fit_mcmc2 <- mcmc(lik_custom7, inits_mcmc, nsteps=100, w=w, upper = upper_mcmc, lower = lower_mcmc)
## 
## # tune w for each parameter
## w <- diff(sapply(fit_mcmc2[2:(ncol(fit_mcmc2)-1)], quantile, c(.05, .95)))
## 
## # run third mcmc for 1000 iter
## fit_mcmc3 <- mcmc(lik_custom7, inits_mcmc, nsteps=1000, w=w, upper = upper_mcmc, lower = lower_mcmc)
## 


## ----plot_mcmc------------------------------------------------------------------------------------------------------------------
# make data long format
d_mcmc <- pivot_longer(fit_mcmc3, names_to = 'param', values_to = 'transition_rate', cols = starts_with('q')) %>%
  left_join(., select(diversitree_df, param, state_1, state_2)) %>%
  left_join(., select(coding, state_1 = hab_pref, state_1_num = hab_pref_num2, state_1_label = initials)) %>%
  left_join(., select(coding, state_2 = hab_pref, state_2_num = hab_pref_num2, state_2_label = initials)) %>%
  mutate(parameter = paste(state_1_label, '->', state_2_label))

# find 95% CIs and bind with ML estimates
d_mcmc_summary <- d_mcmc %>%
  group_by(parameter, state_1, state_2) %>%
  tidybayes::mean_qi(transition_rate) %>%
  left_join(., select(diversitree_df, state_1, state_2, ml_estimate = transition_rate))

# plot ridge plot
ggplot(d_mcmc, aes(transition_rate, forcats::fct_reorder(parameter, transition_rate))) +
  geom_density_ridges() +
  theme_bw(base_size = 14) +
  labs(x = 'transition rate',
       y = 'transition')




## -------------------------------------------------------------------------------------------------------------------------------
# plot 95% CIs
ggplot(d_mcmc_summary, aes(transition_rate, forcats::fct_reorder(parameter, transition_rate))) +
  geom_linerange(aes(xmin = .lower, xmax = .upper)) +
  geom_point(size = 3) +
  geom_point(aes(x = ml_estimate), size = 3, col = 'red') +
  theme_bw(base_size = 14) +
  labs(x = 'transition rate',
       y = 'transition',
       caption = 'red points are ML estimate\nblack points are MCMC average')




## ----mcmc_cor-------------------------------------------------------------------------------------------------------------------
# look at whether the two crazy estimates correlate with each other.
ggplot(fit_mcmc3, aes(q41, q43)) +
  geom_point() +
  theme_bw()



## ----phytools-------------------------------------------------------------------------------------------------------------------

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



## ----simmap---------------------------------------------------------------------------------------------------------------------
## # do stochastic mapping of character traits using phytools
## # feed in the best transition matrix previously found using other methods
## simmap_best <- make.simmap(tree, hab_pref, nsim = 1000, Q = best_matrix)
## 
## # need to split this result up so that files are less than 50MB for GitHub
## 
## # number of splits
## n_splits <- 10
## # find start of each split
## splits <- seq(from = 1, to = 1000, by = 1000/10)
## 
## # save files out
## for(i in 1:length(splits)){
##   # save out every 100 sims
##   saveRDS(simmap_best[splits[i]:(splits[i]+99)], here(paste("data/sequencing_rpoB/processed/transition_rates/simmap/simmap_", i, '.rds', sep = '')))
## }
## 


## ----load_simmap----------------------------------------------------------------------------------------------------------------
# reload simmap files in
simmap_files <- list.files(here("data/sequencing_rpoB/processed/transition_rates/simmap"), full.names = TRUE)
simmap_files <- simmap_files[simmap_files != here('data/sequencing_rpoB/processed/transition_rates/simmap/simmap_summary.rds')]

simmap_best <- purrr::map(simmap_files, readRDS)

simmap_best <- do.call(c, simmap_best)

# load in summary as well!
simmap_summary <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/simmap/simmap_summary.rds'))


## ----summarise_phytools---------------------------------------------------------------------------------------------------------
## # summarise number of switches between states and time spent in each state
## simmap_summary <- describe.simmap(simmap_best, plot=FALSE)
## 
## # remove the tree element - this is the simmap_best
## simmap_summary$tree <- NULL


## ----simmap_summary-------------------------------------------------------------------------------------------------------------
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


## ----common_transitions---------------------------------------------------------------------------------------------------------
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
save_as_image(table_flex, here('plots/sequencing_rpoB/analyses/proportion_of_transitions.png'), zoom = 3, webshot = 'webshot2')

table_flex


## ----plot_best_model------------------------------------------------------------------------------------------------------------
# coerce time into dataframe
d_time <- as.data.frame(simmap_summary$times, col.names = colnames(simmap_summary$times)) %>%
  mutate(n = 1:n()) %>%
  select(-total) %>%
  pivot_longer(cols = c(everything(), -n), names_to = 'state', values_to = 'time') %>%
  group_by(n) %>%
  mutate(prop = time/sum(time)) %>%
  ungroup()

# set colours
cols_hab <- met.brewer('Austria', n = 7)
names(cols_hab) <- c('mud_and_shore', 'freshwater', 'terrestrial', 'freshwater:terrestrial', 'generalist', 'mud_and_shore:terrestrial', 'freshwater:mud_and_shore')
hab_labels <- c('marine mud', 'freshwater', 'terrestrial', 'freshwater + terrestrial', 'generalist', 'marine mud + terrestrial', 'freshwater + marine mud')

# calculate mean time spent in each state
d_timespent <- group_by(d_time, state) %>%
  summarise(mean = mean(prop), .groups = 'drop')

# turn transition matrix into network to plot
d_network <- as_tbl_graph(select(diversitree_df, state_1, state_2, transition_rate)) %>%
  activate(edges) %>%
  filter(!is.na(transition_rate) & transition_rate > 0) %>%
  activate(nodes) %>%
  left_join(., select(d_timespent, name = state, mean)) %>%
  mutate(label = gsub(':', '/', name),
  label = gsub('_', ' ', label),
  label = gsub('mud and shore', 'marine mud', label),
  order = c(1, 2, 6, 7, 3, 4, 5),
  hab1 = gsub(':.*.', '', name),
  hab2 = gsub('.*:', '', name)) %>%
  arrange(order)

p <- ggraph(d_network, layout = 'linear', circular = TRUE) + 
  geom_edge_fan(aes(alpha = transition_rate, 
                width = transition_rate),
                arrow = arrow(length = unit(4, 'mm')),
                end_cap = circle(10, 'mm'),
                start_cap = circle(10, 'mm')) + 
  geom_node_point(aes(size = mean,
                col = name)) +
  theme_void() +
  #geom_node_label(aes(label = label, x=xmin), repel = TRUE) +
  scale_size(range = c(2,20)) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_color_manual('Habitat preference', values = cols_hab, labels = c('terrestrial', 'freshwater', 'marine mud', 'generalist')) +
scale_fill_manual('Habitat preference', values = cols_hab, labels = c('terrestrial', 'freshwater', 'marine mud', 'generalist'))

# grab data for points
point_data <- p$data %>%
  select(x, y, label) %>%
  mutate(nudge_x = ifelse(x < 0, -0.5, 0.5),
  nudge_y = ifelse(y < 0, -0.3, 0.3),
  label = gsub('/', '/\n', label))

p + geom_label(aes(nudge_x + x, nudge_y+y, label = label), point_data, size = MicrobioUoE::pts(18)) +
  theme(legend.position = 'none',
  panel.background = element_rect(fill = 'white', colour = 'white')) +
  coord_cartesian(clip = "off") +
  xlim(c(min(point_data$x) + min(point_data$nudge_x) - 0.3), max(point_data$x) + max(point_data$nudge_x) + 0.3) +
  ylim(c(min(point_data$y) + min(point_data$nudge_y) - 0.1), max(point_data$y) + max(point_data$nudge_y) + 0.1)

# save out model
ggsave(here('plots/sequencing_rpoB/analyses/transition_plot_diversitree.png'), last_plot(), height = 6, width = 8)



## -------------------------------------------------------------------------------------------------------------------------------
# first with the transition rates
d_source_sink_rate <- select(diversitree_df, away = state_1, into = state_2, transition_rate) %>%
  pivot_longer(cols = c(away, into), names_to = 'direction', values_to = 'habitat_preference') %>%
  group_by(habitat_preference, direction) %>%
  summarise(total_rate = sum(transition_rate), .groups = 'drop') %>%
  pivot_wider(names_from = direction, values_from = total_rate) %>%
  mutate(source_sink1 = away / into,
         source_sink2 = into - away)

table_rate <- select(d_source_sink_rate, habitat_preference, away, into, source_sink1) %>%
  mutate(across(away:source_sink1, round, 2),
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
save_as_image(table_rate, here('plots/sequencing_rpoB/analyses/source_sink_rate.png'), zoom = 3, webshot = 'webshot2')

table_rate


## -------------------------------------------------------------------------------------------------------------------------------

# next with simulated counts of transitions on the tree
d_source_sink_count <- select(d_transitions_summary, away = state_1, into = state_2, ave_num) %>%
  pivot_longer(cols = c(away, into), names_to = 'direction', values_to = 'habitat_preference') %>%
  group_by(habitat_preference, direction) %>%
  summarise(total_count = sum(ave_num), .groups = 'drop') %>%
  pivot_wider(names_from = direction, values_from = total_count) %>%
  mutate(source_sink1 = away / into,
         source_sink2 = into - away)
  
table_count <- select(d_source_sink_count, habitat_preference, away, into, source_sink1) %>%
  mutate(across(away:source_sink1, round, 2),
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
save_as_image(table_count, here('plots/sequencing_rpoB/analyses/source_sink_count.png'), zoom = 3, webshot = 'webshot2')

table_count


## ----mass_action_rate-----------------------------------------------------------------------------------------------------------
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
                       across(ave_num:ratio, round, 2),
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
save_as_image(table_expect, here('plots/sequencing_rpoB/analyses/expected_transitions.png'), zoom = 3, webshot = 'webshot2')

table_expect



## ----div_model_prep-------------------------------------------------------------------------------------------------------------

# read in diversification classification
div <- readRDS(here('data/sequencing_rpoB/processed/div_rate_bins.rds'))

# combine d_meta with the diversification rate df
d_meta <- left_join(d_meta, div)

# make sure order of habitat preference is the same in the trait vector as the tip labels
sum(d_meta$tip_label == tree$tip.label) == length(tree$tip.label)
# SUCCESS if TRUE

# create new habitat preference vector that is just saline/non-saline
d_meta <- mutate(d_meta, hab_pref_new = ifelse(str_detect(habitat_preference, 'generalist|mud_and_shore'), 'saline', 'non-saline'),
                 new_trait = paste(hab_pref_new, div_rate, sep = '_'))

# check it is right
group_by(d_meta, new_trait) %>%
  tally()

# replace NA of outgroup - make it freshwater:terrestrial - the most common state
# phytools::make.simmap cannot take NAs
hab_pref_new <- setNames(d_meta$new_trait, d_meta$tip_label)
hab_pref_new2 <- as.numeric(as.factor(hab_pref_new))
names(hab_pref_new2) <- d_meta$tip_label

# create dataframe to easily convert back to actual values
coding_new <- tibble(hab_pref = unname(hab_pref_new), hab_pref_num = unname(hab_pref_new2)) %>%
  distinct() %>%
  arrange(hab_pref) %>%
  mutate(hab_pref_axis = gsub('_', ' ', hab_pref))

coding_new

# check the renaming has worked!
sum(names(hab_pref_new2) == tree$tip.label) == length(tree$tip.label)



## ----saline_diversitree_noshow--------------------------------------------------------------------------------------------------
# read in models
mod_er_div <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/mod_er_div.rds'))
mod_ard_div <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/mod_ard_div.rds'))
mod_sym_div <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/mod_sym_div.rds'))

# all of these methods needs a likelihood function, we can build a Mkn model
lik_ard_div <- make.mkn(tree, hab_pref_new2, k = max(hab_pref_new2))


## ----div_diversitree_norun------------------------------------------------------------------------------------------------------
## # all of these methods needs a likelihood function, we can build a Mkn model
## lik_ard_div <- make.mkn(tree, hab_pref_new2, k = max(hab_pref_new2))
## 
## argnames(lik_ard_div)
## 
## # make symmetric rates model
## lik_sym_div <- constrain(lik_ard_div,
##                      q12~q21, q13~q31, q14~q41,
##                      q23~q32, q24~q42,
##                      q34~q43
##                      )
## 
## argnames(lik_sym_div)
## 
## # make equal rates model
## lik_er_div <- constrain(lik_ard_div,
##                     q13~q12, q14~q12,
##                     q21~q12, q23~q12, q24~q12,
##                     q31~q12, q32~q12, q34~q12,
##                     q41~q12, q42~q12, q43~q12)
## 
## argnames(lik_er_div)
## 
## # need to pass start values to it - can grab these from the ape::ace, but we will just pass an average rate to the model
## inits_ard <- rep(1, length(argnames(lik_ard_div)))
## inits_sym <- rep(1, length(argnames(lik_sym_div)))
## inits_er <- rep(1, length(argnames(lik_er_div)))
## 
## # find the maximum likelihood estimates of this model
## mod_er_div <- find.mle(lik_er_div, inits_er, method = 'subplex', control = list(maxit = 50000))
## mod_sym_div <- find.mle(lik_sym_div, inits_sym, method = 'subplex', control = list(maxit = 50000))
## mod_ard_div <- find.mle(lik_ard_div, inits_ard, method = 'subplex', control = list(maxit = 50000))


## -------------------------------------------------------------------------------------------------------------------------------
# check AIC
AIC(mod_er_div, mod_sym_div, mod_ard_div) %>% 
  arrange(AIC)



## ----plot_div-------------------------------------------------------------------------------------------------------------------
get_diversitree_df(mod_ard_div, coding_new$hab_pref_num, coding_new$hab_pref) %>%
  left_join(., select(coding_new, state_1 = hab_pref, state_1_num = hab_pref_num, state_1_label = hab_pref_axis)) %>%
  left_join(., select(coding_new, state_2 = hab_pref, state_2_num = hab_pref_num, state_2_label = hab_pref_axis)) %>%
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
  labs(y = 'state 1',
  x = 'state 2',
  title = paste('all rates different with', length(mod_ard_div$par), 'free parameters', sep = ' ')) +
  coord_fixed() +
  scale_color_manual(values = c('red', 'black'))

ggsave(here('plots/sequencing_rpoB/analyses/transition_matrix_fast_slow.png'), last_plot(), height = 6, width = 8)



## ----setup_new_trait------------------------------------------------------------------------------------------------------------

# create new habitat preference vector that is just saline/non-saline generalist / non generalist
d_meta <- mutate(d_meta, saline_or_no = ifelse(str_detect(habitat_preference, 'generalist|mud_and_shore'), 'saline', 'non-saline'),
                 specialist_or_no = ifelse(str_detect(habitat_preference, ':|generalist'), 'generalist', 'specialist'),
                 new_trait = paste(saline_or_no, specialist_or_no, sep = '_'))

# check it is right
group_by(d_meta, new_trait) %>%
  tally()

# replace NA of outgroup - make it freshwater:terrestrial - the most common state
# phytools::make.simmap cannot take NAs
hab_pref_new <- setNames(d_meta$new_trait, d_meta$tip_label)
hab_pref_new2 <- as.numeric(as.factor(hab_pref_new))
names(hab_pref_new2) <- d_meta$tip_label

# create dataframe to easily convert back to actual values
coding_new <- tibble(hab_pref = unname(hab_pref_new), hab_pref_num = unname(hab_pref_new2)) %>%
  distinct() %>%
  arrange(hab_pref) %>%
  mutate(hab_pref_axis = gsub('_', ' ', hab_pref))

coding_new

# check the renaming has worked!
sum(names(hab_pref_new2) == tree$tip.label) == length(tree$tip.label)


## ----premusse_ard_noshow--------------------------------------------------------------------------------------------------------
# read in models
mod_ard_premusse <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/mod_ard_premusse.rds'))
mod_sym_premusse <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/mod_sym_premusse.rds'))
mod_premusse1 <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/mod_premusse1.rds'))
mod_premusse2 <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/mod_premusse2.rds'))

# set up likelihood model for diversitree
lik_premusse <- make.mkn(tree, hab_pref_new2, k = max(hab_pref_new2))

# make custom matrix model
lik_premusse1 <- constrain(lik_premusse, 
                   q14~0, q23~0, q41~0)


## ----do_ard_selection-----------------------------------------------------------------------------------------------------------
## # set up likelihood model for diversitree
## lik_premusse <- make.mkn(tree, hab_pref_new2, k = max(hab_pref_new2))
## 
## lik_sym_premusse <- constrain(lik_premusse,
##                      q12 ~ q21, q13 ~ q31, q14 ~ q41,
##                      q23 ~ q32, q24 ~ q42,
##                      q34 ~ q43)
## 
## inits <- rep(1, length(argnames(lik_premusse)))
## inits_sym <- rep(1, length(argnames(lik_sym_premusse)))
## 
## argnames(lik_sym)
## 
## # run ARD model
## mod_ard_premusse <- find.mle(lik, inits_ard, method = 'subplex', control = list(maxit = 50000))
## mod_sym_premusse <- find.mle(lik_sym, inits_sym, method = 'subplex', control = list(maxit = 50000))
## 


## ----do_ard_selection2----------------------------------------------------------------------------------------------------------
# check AIC
AIC(mod_ard_premusse, mod_sym_premusse) %>%
  arrange(AIC)

# look at estimates
mod_ard_premusse$par %>% sort()

# make custom matrix model
lik_premusse1 <- constrain(lik_premusse, 
                   q14~0, q23~0, q41~0)

# make start parameters
inits_custom1 <- rep(1, length(argnames(lik_premusse1)))


## -------------------------------------------------------------------------------------------------------------------------------
## # run model
## mod_premusse1 <- find.mle(lik_premusse1, inits_custom1, method = 'subplex', control = list(maxit = 50000))


## -------------------------------------------------------------------------------------------------------------------------------
# do model selection
AIC(mod_ard_premusse, mod_premusse1) %>%
  arrange(AIC)

anova(mod_ard_premusse, mod_premusse1)

mod_premusse1$par %>% sort()

# make custom matrix model
lik_premusse2 <- constrain(lik_premusse1, 
                   q24~0)

# make start parameters
inits_custom2 <- rep(1, length(argnames(lik_premusse2)))


## -------------------------------------------------------------------------------------------------------------------------------
## # run model
## mod_premusse2 <- find.mle(lik_2, inits_custom2, method = 'subplex', control = list(maxit = 50000))
## 


## -------------------------------------------------------------------------------------------------------------------------------

# do ANOVAs
anova(mod_premusse1, mod_premusse2)

# check AIC
AIC(mod_ard_premusse, mod_premusse1, mod_premusse2) %>%
  arrange(AIC)


## ----plot_matrix_again----------------------------------------------------------------------------------------------------------
# code to plot transition matrix
get_diversitree_df(mod_premusse1, coding_new$hab_pref_num, coding_new$hab_pref) %>%
  left_join(., select(coding_new, state_1 = hab_pref, state_1_num = hab_pref_num, state_1_label = hab_pref_axis)) %>%
  left_join(., select(coding_new, state_2 = hab_pref, state_2_num = hab_pref_num, state_2_label = hab_pref_axis)) %>%
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
  labs(y = 'state 1',
  x = 'state 2',
  title = paste('all rates different with', length(mod_premusse1$par), 'free parameters', sep = ' ')) +
  coord_fixed() +
  scale_color_manual(values = c('red', 'black'))

ggsave(here('plots/sequencing_rpoB/analyses/transition_matrix_premusse.png'), last_plot(), height = 6, width = 8)



## ----fit_musse_noshow-----------------------------------------------------------------------------------------------------------
# read in models
fit_musse <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/mod_musse_full.rds'))
fit_musse_no_se <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/mod_musse_no_se.rds'))
fit_musse_no_ss <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/mod_musse_no_ss.rds'))
fit_musse_no_se2 <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/mod_musse_no_se2.rds'))
fit_musse_no_sse <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/mod_musse_no_sse.rds'))
fit_mcmc3 <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/fit_mcmc_musse.rds'))



## ----fit_musse------------------------------------------------------------------------------------------------------------------
## # set up sampling fractions, set them all to 1
## sampling_frac <- setNames(rep(1, times = 4), sort(unique(hab_pref_new2)))
## 
## # set up likelihood model for diversitree
## lik_musse <- make.musse(tree, hab_pref_new2, k = max(hab_pref_new2), sampling.f = sampling_frac)
## 
## # set constraints for transitions that do not occur
## lik_musse <- constrain(lik_musse, q14~0, q23~0, q41~0)
## 
## # we can estimate starting values using starting.point.musse()
## start_vals <- starting.point.musse(tree, k = max(hab_pref_new2))
## 
## for(i in 1:length(mod_1$par)){
##     start_vals[names(start_vals) == names(mod_premusse1$par)[i]] <- unname(mod_premusse1$par[i])
## }
## 
## # fit musse model
## fit_musse <- find.mle(lik_musse, x.init = start_vals[argnames(lik_musse)], method = 'subplex', control = list(maxit = 50000))
## 


## ----run_null_models------------------------------------------------------------------------------------------------------------
## # remove state dependent speciation and extinction parameters
## lik_null_no_sse <- constrain(lik_musse, lambda2 ~ lambda1, lambda3 ~ lambda1, lambda4 ~ lambda1,
##                              mu2 ~ mu1, mu3 ~ mu1, mu4 ~ mu1)
## 
## # remove only speciation parameters
## lik_null_no_ss <- constrain(lik_musse, lambda2 ~ lambda1, lambda3 ~ lambda1, lambda4 ~ lambda1)
## 
## # remove only extinction parameters
## lik_null_no_se <- constrain(lik_musse, mu2 ~ mu1, mu3 ~ mu1, mu4 ~ mu1)
## 
## # fit model
## fit_musse_no_sse <- find.mle(lik_null_no_sse, x.init = start_vals[argnames(lik_null_no_sse)], method = 'subplex', control = list(maxit = 50000))
## fit_musse_no_ss <- find.mle(lik_null_no_ss, x.init = start_vals[argnames(lik_null_no_ss)],method = 'subplex', control = list(maxit = 50000))
## fit_musse_no_se <- find.mle(lik_null_no_se, x.init = start_vals[argnames(lik_null_no_se)], method = 'subplex', control = list(maxit = 50000))
## 


## ----musse_anova----------------------------------------------------------------------------------------------------------------
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


## ----rerun_best-----------------------------------------------------------------------------------------------------------------
## # we can estimate starting values using starting.point.musse()
## start_vals <- starting.point.musse(tree, k = max(hab_pref_new2))
## 
## fit_musse_no_se2 <- find.mle(lik_null_no_se, x.init = start_vals[argnames(lik_null_no_se)], method = 'subplex', control = list(maxit = 50000))


## ----compare_sse_models---------------------------------------------------------------------------------------------------------

AIC(fit_musse_no_se, fit_musse_no_se2)

fit_musse_no_se$lnLik
fit_musse_no_se2$lnLik



## -------------------------------------------------------------------------------------------------------------------------------

fit_musse_no_se$par %>% round(2)

coding_new


## ----musse_mcmc-----------------------------------------------------------------------------------------------------------------
## # set up initial start values
## inits_mcmc <- fit_musse_no_se$par
## 
## # set up upper and lower limits - limit the values to be < 10 times the max value
## lower_mcmc <- rep(0, length(inits_mcmc))
## upper_mcmc <- rep(max(inits_mcmc)*100, length(inits_mcmc))
## 
## # run first mcmc to tune w
## fit_mcmc <- mcmc(lik_null_no_se, inits_mcmc, nsteps = 10, w = 0.1, upper = upper_mcmc, lower = lower_mcmc)
## 
## # tune w for each parameter
## w <- diff(sapply(fit_mcmc[2:(ncol(fit_mcmc)-1)], quantile, c(.05, .95)))
## 
## # run second mcmc to tune w
## fit_mcmc2 <- mcmc(lik_null_no_se, inits_mcmc, nsteps=100, w=w, upper = upper_mcmc, lower = lower_mcmc)
## 
## # tune w for each parameter
## w <- diff(sapply(fit_mcmc2[2:(ncol(fit_mcmc2)-1)], quantile, c(.05, .95)))
## 
## # run third mcmc for 1000 iter
## fit_mcmc3 <- mcmc(lik_null_no_se, inits_mcmc, nsteps=1000, w=w, upper = upper_mcmc, lower = lower_mcmc)
## 


## ----plot_speciation------------------------------------------------------------------------------------------------------------
# set colours
cols_hab <- met.brewer('Ingres', n = 4)
names(cols_hab) <- c('saline_specialist', 'saline_generalist', 'non-saline_generalist', 'non-saline_specialist')
hab_labels <- c('saline specialist', 'saline generalist', 'non-saline generalist', 'non-saline specialist')

# first grab the mcmc results for speciation
d_speciation <- select(fit_mcmc3, i, starts_with('lambda')) %>%
  pivot_longer(cols = starts_with('lambda'), names_to = 'hab_pref_num', values_to = 'speciation_rate', names_prefix = 'lambda') %>%
  mutate(hab_pref_num = as.numeric(hab_pref_num)) %>%
  left_join(., coding_new) %>%
  separate(hab_pref, c('habitat', 'type'), sep = '_', remove = FALSE)

d_speciation_summary <- group_by(d_speciation, hab_pref_num, hab_pref, habitat, type, hab_pref_axis) %>%
  tidybayes::mean_qi(speciation_rate)

# grab out the maximum likelihood estimates
d_speciation_ml <- fit_musse_no_se$par[grepl('lambda', names(fit_musse_no_se$par))] %>%
  data.frame(speciation_rate = ., hab_pref_num = parse_number(names(.)), row.names = NULL) %>%
  left_join(., coding_new) %>%
  separate(hab_pref, c('habitat', 'type'), sep = '_', remove = FALSE)

# plot
ggplot(d_speciation_summary, aes(hab_pref_axis, speciation_rate, col = hab_pref)) +
  geom_point(size = 8, show.legend = FALSE) +
  geom_linerange(aes(ymin = .lower, ymax = .upper), show.legend = FALSE) +
  geom_point(aes(hab_pref_axis, speciation_rate), col = 'grey', data = d_speciation_ml) +
  theme_bw(base_size = 16) +
  scale_color_manual(values = cols_hab) +
  ylim(c(0,16)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  labs(x = 'Habitat preference',
       y = 'Speciation rate')
  
ggsave(here('plots/sequencing_rpoB/analyses/musse_speciation.png'), last_plot(), height = 5, width = 7)
  


## ----plot_musse_transitions-----------------------------------------------------------------------------------------------------

# grab out transition rates
d_transition_musse <- select(fit_mcmc3, i, starts_with('q')) %>%
  pivot_longer(cols = starts_with('q'), names_to = 'transition', values_to = 'transition_rate', names_prefix = 'q') %>%
  mutate(state_1_num = as.numeric(substr(transition, 1,1)),
         state_2_num = as.numeric(substr(transition, 2,2))) %>%
  left_join(., select(coding_new, state_1_num = hab_pref_num, state_1 = hab_pref)) %>%
  left_join(., select(coding_new, state_2_num = hab_pref_num, state_2 = hab_pref)) %>%
  group_by(state_1_num, state_1, state_2_num, state_2) %>%
  tidybayes::mean_qi(transition_rate)

# check it is right
d_size <- group_by(d_meta, new_trait) %>%
  tally() %>%
  mutate(prop = n/sum(n))

# turn transition matrix into network to plot
d_network <- as_tbl_graph(select(d_transition_musse, state_1, state_2, transition_rate)) %>%
  activate(edges) %>%
  filter(!is.na(transition_rate) & transition_rate > 0) %>%
  activate(nodes) %>%
  left_join(., select(d_size, name = new_trait, prop)) %>%
  mutate(label = name,
         label = gsub('_', ' ', label),
  order = c(1, 2, 3, 4)) %>%
  arrange(order)

p <- ggraph(d_network, layout = 'linear', circular = TRUE) + 
  geom_edge_fan(aes(alpha = transition_rate, 
                width = transition_rate),
                arrow = arrow(length = unit(4, 'mm')),
                end_cap = circle(10, 'mm'),
                start_cap = circle(10, 'mm'),
                show.legend = FALSE) + 
  geom_node_point(aes(size = prop,
                col = name),
                show.legend = FALSE) +
  theme_void() +
  scale_size(range = c(2,20)) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_color_manual(values = cols_hab)

# grab data for points
point_data <- p$data %>%
  select(x, y, label) %>%
  mutate(nudge_x = ifelse(x < 0, -0.5, 0.5),
  nudge_y = ifelse(y < 0, -0.3, 0.3),
  label = gsub(' ', '\n', label))

p + geom_label(aes(nudge_x + x, nudge_y+y, label = label), point_data, size = MicrobioUoE::pts(18)) +
  theme(legend.position = 'none',
  panel.background = element_rect(fill = 'white', colour = 'white')) +
  coord_cartesian(clip = "off") +
  xlim(c(min(point_data$x) + min(point_data$nudge_x) - 0.3), max(point_data$x) + max(point_data$nudge_x) + 0.3) +
  ylim(c(min(point_data$y) + min(point_data$nudge_y) - 0.1), max(point_data$y) + max(point_data$nudge_y) + 0.1)

# save out model
ggsave(here('plots/sequencing_rpoB/analyses/transition_plot_musse.png'), last_plot(), height = 6, width = 8)


## ----setup muhisse--------------------------------------------------------------------------------------------------------------
# make new coding
hab_pref_new3 <- case_when(hab_pref_new2 == 1 ~ '00',
                           hab_pref_new2 == 2 ~ '01',
                           hab_pref_new2 == 3 ~ '10',
                           hab_pref_new2 == 4 ~ '11')
names(hab_pref_new3) <- names(hab_pref_new2)

coding_new <- mutate(coding_new, muhisse_coding = c('00', '01', '10', '11'))

# check the renaming has worked!
sum(names(hab_pref_new3) == tree$tip.label) == length(tree$tip.label)

# create transition rate matrix
trans_mat <- TransMatMakerMuHiSSE(hidden.traits = 1)

# we now need to customise some of the transition rates and make some of them possible that are not currently
# namely 10 (saline generalism) to 01 (non-saline specialist)
trans_mat['(10A)', '(01A)'] <- max(trans_mat, na.rm = TRUE) + 1
trans_mat['(10B)', '(01B)'] <- max(trans_mat, na.rm = TRUE) + 1

# As the best MuSSE model just has variation in speciation, this is what we will fit here
turnover <- c(1,2,3,4,5,6,7,8)
extinction_fraction <- rep(1, 8) 

# set estimated proportion of extant species
f = c(1,1,1,1)

# get the dataset ready
muhisse_df <- data.frame(tip_label = names(hab_pref_new3),
                         saline = substr(hab_pref_new3, 1,1),
                         specialist = substr(hab_pref_new3, 2,2))


## ----run_MuHiSSE----------------------------------------------------------------------------------------------------------------
## muhisse_1 <- MuHiSSE(phy=tree, data=muhisse_df, f=f,
##                      turnover=turnover,
##                      eps=extinction_fraction,
##                      hidden.states=TRUE,
##                      trans.rate=trans_mat)


## ----old_code-------------------------------------------------------------------------------------------------------------------
## # first get the trait matrix into a data frame we can use
## trans_mat_df <- as_tibble(trans_mat) %>%
##   mutate(state1_full = row.names(trans_mat)) %>%
##   pivot_longer(cols = starts_with('('), values_to = 'param_num', names_to = 'state2_full') %>%
##   mutate(state1 = substr(state1_full, 2,3),
##          state2 = substr(state2_full, 2,3)) %>%
##   left_join(., get_diversitree_df(mod_premusse1, trait_vec = coding_new$hab_pref_num, replace_vec = coding_new$muhisse_coding))
## 
## premusse_rates <- get_diversitree_df(mod_premusse1, trait_vec = coding_new$hab_pref_num, replace_vec = coding_new$muhisse_coding) %>%
##   select(state1 = state_1, state2 = state_2, transition_rate)
## 
## trans_mat_df <- left_join(trans_mat_df, premusse_rates)
## 

