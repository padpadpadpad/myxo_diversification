
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

mod_sym <- find.mle(lik_sym, inits_sym, method = 'subplex', control = list(maxit = 50000))
mod_er <- find.mle(lik_er, inits_er, method = 'subplex', control = list(maxit = 50000))
mod_ard <- find.mle(lik_ard, inits_ard, method = 'subplex', control = list(maxit = 50000))

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

# make custom matrix model
lik_custom6 <- constrain(lik_ard, 
                     q13~0, q15~0, q16~0, q17~0, q26~0, q45~0, q51~0, q53~0, q54~0, q62~0,
                     q31~0, q35~0, q36~0, q63~0, q72~0, q73~0, q75~0,
                     q57~0,
                     q43~0,
                     q56~0,
                     q76~0)

# make start parameters
inits_custom6 <- rep(1, length(argnames(lik_custom6)))

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

# make custom matrix model
lik_custom7 <- constrain(lik_ard, 
                         q13~0, q15~0, q16~0, q17~0, q26~0, q45~0, q51~0, q53~0, q54~0, q62~0,
                         q31~0, q35~0, q36~0, q63~0, q72~0, q73~0, q75~0,
                         q57~0,
                         q43~0,
                         q56~0,
                         q76~0,
                         q64~0)

# make start parameters
inits_custom7 <- rep(1, length(argnames(lik_custom7)))

## # run model
mod_custom7 <- find.mle(lik_custom7, inits_custom7, method = 'subplex', control = list(maxit = 50000))


## ----diversitree_7_compare------------------------------------------------------------------------------------------------------
# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_custom1, mod_custom2, mod_custom3, mod_custom4, mod_custom5, mod_custom6, mod_custom7) %>% arrange(AIC)

# anova
anova(mod_custom6, mod_custom7)

# sort parameter estimates
mod_custom7$par %>% sort()

## ----diversitree_8_load--------------------------------------------------------------------------------------

# make custom matrix model
lik_custom8 <- constrain(lik_ard, 
                        q13~0, q15~0, q16~0, q17~0, q26~0, q45~0, q51~0, q53~0, q54~0, q62~0,
                        q31~0, q35~0, q36~0, q63~0, q72~0, q73~0, q75~0,
                        q57~0,
                        q43~0,
                        q56~0,
                        q76~0,
                        q64~0,
                        q65~0)

# make start parameters
inits_custom8 <- rep(1, length(argnames(lik_custom8)))

# run model
mod_custom8 <- find.mle(lik_custom8, inits_custom8, method = 'subplex', control = list(maxit = 50000))


## ----diversitree_8_compare------------------------------------------------------------------------------------------------------
# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_custom1, mod_custom2, mod_custom3, mod_custom4, mod_custom5, mod_custom6, mod_custom7, mod_custom8) %>% arrange(AIC)

# anova
anova(mod_custom7, mod_custom8)

# sort parameter estimates
mod_custom8$par %>% sort()

# make custom matrix model
lik_custom9 <- constrain(lik_ard, 
                         q13~0, q15~0, q16~0, q17~0, q26~0, q45~0, q51~0, q53~0, q54~0, q62~0,
                         q31~0, q35~0, q36~0, q63~0, q72~0, q73~0, q75~0,
                         q57~0,
                         q43~0,
                         q56~0,
                         q76~0,
                         q64~0,
                         q65~0,
                         q46~0)

# make start parameters
inits_custom9 <- rep(1, length(argnames(lik_custom9)))

# run model
mod_custom9 <- find.mle(lik_custom9, inits_custom9, method = 'subplex', control = list(maxit = 50000))


## ----diversitree_9_compare------------------------------------------------------------------------------------------------------
# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_custom1, mod_custom2, mod_custom3, mod_custom4, mod_custom5, mod_custom6, mod_custom7, mod_custom8, mod_custom9) %>% arrange(AIC)

# anova
anova(mod_custom8, mod_custom9)

# sort parameter estimates
mod_custom9$par %>% sort()

# mod custom 4 is the best.
# run it for more iterations
mod_custom4 <- find.mle(lik_custom4, inits_custom4, method = 'subplex', control = list(maxit = 100000))

# save out models so far
saveRDS(mod_ard, 'data/sequencing_rpoB/processed/transition_rates/asv_mod_ard.rds')
saveRDS(mod_sym, 'data/sequencing_rpoB/processed/transition_rates/asv_mod_sym.rds')
saveRDS(mod_er, 'data/sequencing_rpoB/processed/transition_rates/asv_mod_er.rds')
saveRDS(mod_custom1, 'data/sequencing_rpoB/processed/transition_rates/asv_mod_custom1.rds')
saveRDS(mod_custom2, 'data/sequencing_rpoB/processed/transition_rates/asv_mod_custom2.rds')
saveRDS(mod_custom3, 'data/sequencing_rpoB/processed/transition_rates/asv_mod_custom3.rds')
saveRDS(mod_custom4, 'data/sequencing_rpoB/processed/transition_rates/asv_mod_custom4.rds')
saveRDS(mod_custom5, 'data/sequencing_rpoB/processed/transition_rates/asv_mod_custom5.rds')
saveRDS(mod_custom6, 'data/sequencing_rpoB/processed/transition_rates/asv_mod_custom6.rds')
saveRDS(mod_custom7, 'data/sequencing_rpoB/processed/transition_rates/asv_mod_custom7.rds')
saveRDS(mod_custom8, 'data/sequencing_rpoB/processed/transition_rates/asv_mod_custom8.rds')
saveRDS(mod_custom9, 'data/sequencing_rpoB/processed/transition_rates/asv_mod_custom9.rds')


## ----plot_transition_matrix-----------------------------------------------------------------------------------------------------
diversitree_df <- get_diversitree_df(mod_custom4, coding$hab_pref_num2, coding$hab_pref)

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
  title = paste('all rates different with', length(mod_custom4$par), 'free parameters', sep = ' ')) +
  coord_fixed() +
  scale_color_manual(values = c('red', 'black'))


## ----mcmc_try-------------------------------------------------------------------------------------------------------------------
## # set up initial start values
inits_mcmc <- mod_custom4$par

## # set up upper and lower limits - limit the values to be < 10 times the max value
lower_mcmc <- rep(0, length(inits_mcmc))
upper_mcmc <- rep(max(inits_mcmc)*10, length(inits_mcmc))
## 
## # run first mcmc to tune w
fit_mcmc <- mcmc(lik_custom4, inits_mcmc, nsteps = 10, w = 0.1, upper = upper_mcmc, lower = lower_mcmc)

## # tune w for each parameter
w <- diff(sapply(fit_mcmc[2:(ncol(fit_mcmc)-1)], quantile, c(.05, .95)))

# run second mcmc to tune w
fit_mcmc2 <- mcmc(lik_custom4, inits_mcmc, nsteps=100, w=w, upper = upper_mcmc, lower = lower_mcmc)

## # tune w for each parameter
w <- diff(sapply(fit_mcmc2[2:(ncol(fit_mcmc2)-1)], quantile, c(.05, .95)))

# run third mcmc for 1000 iter
fit_mcmc3 <- mcmc(lik_custom4, inits_mcmc, nsteps=1000, w=w, upper = upper_mcmc, lower = lower_mcmc)


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


