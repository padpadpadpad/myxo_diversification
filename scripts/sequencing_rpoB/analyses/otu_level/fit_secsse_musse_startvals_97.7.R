# fit the MuSSE model using secsse

# load in packages ####

# make sure curl is installed
library(curl)
librarian::shelf(diversitree, secsse, DDD, doParallel, foreach, doMC, tidyverse, here, furrr)

# identify conflicts in the tidyverse packages and other packages
tidyverse_conflicts()

# load in data ####

# load in data ####

# filename
name <- 'musse_97.7'

# estimate all transition rates or just the hidden transitions
# to reduce number of estimatable parameters we fixed the transition rates of transitions from the markov model
# can be "hidden" or "all"
transitions <- 'hidden'

# server - yes or no
server <- FALSE

if(server == FALSE){
  # read in habitat preference
  d_habpref <- read.csv(here('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_97.7percent.csv'))
  # read in phyloseq object and grab tax table
  d_taxa <- readRDS(here('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_97.7percent_filt.rds'))
  # read in tree
  tree <-  read.tree('data/sequencing_rpoB/raxml/trees/myxo_97.7/myxo_97.7_treepl_cv_node_labels.tre')
  # read in Mk model
  fit_mk <- readRDS('data/sequencing_rpoB/processed/transition_rates/mod_custom_3_otu97.7.rds')
}

# alter tip labels to remove family as they will not link to the distance matrix
# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

tree$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

d_taxa <- d_taxa %>%
  phyloseq::tax_table() %>%
  data.frame() %>%
  janitor::clean_names() %>%
  rownames_to_column('otu')

# collapse full generalists, marine + terrestrial generalist, and freshwater + marine generalist into a single preference
# create new habitat preference vector
d_habpref <- mutate(d_habpref, habitat_preference2 = ifelse(habitat_preference %in% c('freshwater:marine_mud:terrestrial', 'freshwater:marine_mud', 'marine_mud:terrestrial'), 'marine_mud_generalist', habitat_preference),
                    # rename all habitat preference vectors for easy renaming
                    habitat_preference3 = case_when(habitat_preference2 == 'marine_mud_generalist' ~ "marine mud generalist",
                                                    habitat_preference2 == 'marine_mud' ~ "marine mud specialist",
                                                    habitat_preference == "terrestrial" ~ "terrestrial specialist",
                                                    habitat_preference2 == "freshwater" ~ "freshwater specialist",
                                                    habitat_preference2 == "freshwater:terrestrial" ~ "freshwater + terrestrial generalist"))

# create d_meta
# use habitat_preference3 which collapses into marine mud generalists
d_meta <- left_join(dplyr::select(d_habpref, otu, habitat_preference = habitat_preference3, num_present), dplyr::select(d_taxa, otu:family))

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

# try and run SecSSE which runs concealed state and speciation models
# https://cran.r-project.org/web/packages/secsse/vignettes/Using_secsse.html

# set up secsse model
trait <- data.frame(otu = names(hab_pref_num), trait = hab_pref_num)
traits <- sortingtraits(trait, tree)

# setup arguments to pass to secsse_ml

# set number of concealed states
num_concealed_states <- 1

# setup parameter list
idparslist <- id_paramPos(traits, num_concealed_states = num_concealed_states)

idparslist

# setup speciation rates ####
idparslist$lambdas[] <- c(1,2,3,4,5)

# setup extinction rates ####
# firstly make all extinction rates the same
idparslist$mus[] <- 6

# setup transition rates ####

# make a bunch of transitions 0 
# these transitions were not possible in the markov model
# q14~0, q24~0, q25~0, q41~0, q53~0, q42~0, q45~0, q54~0, q52~0

# make transition matrix a dataframe so I can set rules more easily
q_matrix <- data.frame(idparslist$Q) %>%
  rownames_to_column(var = 'from') %>%
  pivot_longer(starts_with('X'), values_to = 'id', names_to = 'to') %>%
  mutate(to = gsub('X', '', to),
         from_trait = parse_number(from),
         from_hidden = substr(from,2,2),
         to_trait = substr(to, 1,1),
         to_hidden = substr(to, 2,2),
         transition = paste('q', from_trait, to_trait, sep = ''),
         id = -id) 

# impossible transitions
zero_transitions <- fit_mk$par.full[fit_mk$par.full == 0] %>% names()

# first make any of the transitions not possible in the Markov model 0
q_matrix <- mutate(q_matrix,
                   new_id = ifelse(transition %in% zero_transitions, 0, id))

# run a for loop to replace each number in the initial q matrix
q <- idparslist[[3]]

for(i in min(q,na.rm = TRUE):max(q, na.rm = TRUE)){
  q[which(q == i)] <- filter(q_matrix, -id == i) %>% pull(new_id)
}

q

idparslist$Q <- q

# replace all values in q matrix that are not 0 or NA by next logical value

# find the value for the extinction parameter
mu <- unique(idparslist$mus)

# find number of params
q_params <- idparslist$Q[!is.na(idparslist$Q) & idparslist$Q != 0] %>% unique() %>% length()

# create index data frame to change current values to next logical value
to_change <- data.frame(val = idparslist$Q[!is.na(idparslist$Q) & idparslist$Q != 0] %>% unique(),
                        new = (mu+1):(q_params+mu))

for(i in 1:nrow(to_change)){
  idparslist$Q[which(idparslist$Q == to_change$val[i])] <- to_change$new[i]
}

idparslist$Q

# set initial values ####

# mk transition rates
mk_transitions <- fit_mk$par[grepl('q', names(fit_mk$par))]

# for lambda and mu
startingpoint <- bd_ML(brts = ape::branching.times(tree))
init_lambda <- startingpoint$lambda0
init_mu <- startingpoint$mu0

# set initial values for transition rates based on those from diversitree
init_transition <- data.frame(idparslist$Q) %>%
  rownames_to_column(var = 'from') %>%
  pivot_longer(starts_with('X'), values_to = 'id', names_to = 'to') %>%
  mutate(transition = paste('q', parse_number(from), parse_number(to), sep = '')) %>%
  filter(id > 0) %>%
  arrange(id) %>%
  select(id, transition) %>%
  distinct() %>%
  left_join(., data.frame(transition = names(mk_transitions), rate = unname(mk_transitions))) %>%
  mutate(rate = replace_na(rate, mean(rate, na.rm = TRUE))) %>%
  select(id, rate) %>%
  distinct()

initparsopt <- c(rep(init_lambda, times = max(idparslist$lambdas)),
                 rep(init_mu, times = 1))

# check number of estimated parameters is the same as number of initial values
idparsopt <- c(1:max(idparslist$Q, na.rm=TRUE))

idparslist

# set the ID and values for the fixed parameters
# fix the values of the transition rates we estimated from the Mk model
idparsfix <- c(0, init_transition$id) # zeroes have the value of zero
parsfix <- c(0, init_transition$rate)

# remove any parameters from the idparsopt (to estimate) that are present in idparsfix (parameters with fixed values)
idparsopt <- idparsopt[!idparsopt %in% idparsfix]

length(initparsopt) == length(idparsopt)

# set number of iterations
max_iter <- 1000 * round((1.25)^length(idparsopt))

# test initial values to work out their initial log likelihood
idparslist

# write function to get initial values into the correct format
get_inits_matrix <- function(inits, idparsopt, idparsfix, parsfix){
  
  id_fixed <- parsfix
  names(id_fixed) <- idparsfix
  
  for(i in 1:length(idparslist)){
    for(j in 1:length(idparslist[[i]])){
      if(idparslist[[i]][j] %in% 1:length(inits)){idparslist[[i]][j] <- inits[idparslist[[i]][j]]}
      else next 
    }
    for(j in 1:length(idparslist[[i]])){
      if(idparslist[[i]][j] %in% names(id_fixed)){idparslist[[i]][j] <- id_fixed[names(id_fixed) == idparslist[[i]][j]]}
      else next 
    }
  }
  return(idparslist)
}

# set up parameter combinations

# name the start values
inits_lambda <- rep(init_lambda, times = max(idparslist$lambdas))
inits_mu <- rep(init_mu, times = length(unique(idparslist$mus)))
inits_q <- init_transition

# create multiplication factors for the initial values

# 50% and double each set
vals <- c(0.5,1,2)

# create full grid of start values
# add in q here if wanting to estimate the transition rates alongside speciation and extinction rates
vals <- expand.grid(lambda = vals, mu = vals) %>%
  as.tibble()

# create an empty dataframe
inits_ml <- mutate(vals, 
                   inits = list(NA), 
                   loglik = NA,
                   id = 1:n())

# set up for loop to run screen for initial values
pb <- progress::progress_bar$new(total = nrow(inits_ml))

for(i in 1:nrow(inits_ml)){
  pb$tick()
  
  # create inits
  temp_inits <- c((inits_lambda*inits_ml$lambda[i]),
                  (inits_mu*inits_ml$mu[i]))
  
  temp_mat <- get_inits_matrix(temp_inits, idparslist, idparsfix, parsfix)
  
  temp_ml <- # check maximum likelihood values of initial values
    secsse_loglik(
      temp_mat,
      tree,
      traits,
      num_concealed_states = num_concealed_states,
      cond = "maddison_cond",
      root_state_weight = "maddison_weights",
      sampling_fraction = rep(1, times = length(unique(traits))),
      see_ancestral_states = FALSE
    )
  
  inits_ml$inits[[i]] <- temp_inits
  
  inits_ml$loglik[i] <- temp_ml
  
}

# save this out
saveRDS(inits_ml, paste('data/sequencing_rpoB/processed/secsse/init_vals_ml/otu_cutoffs/', name, '.rds', sep = ''))

ggplot(inits_ml, aes(loglik)) +
  geom_histogram()
