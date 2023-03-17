# try a five hidden state concealed trait model using secSSE

# load in packages ####

# make sure curl is installed
library(curl)
librarian::shelf(diversitree, secsse, DDD, apTreeshape, doParallel, foreach, doMC, tidyverse, here)

# identify conflicts in the tidyverse packages and other packages
tidyverse_conflicts()

# load in data ####

# filename
name <- 'twohidden_notraitspeciation'

# server - yes or no
server <- TRUE

if(server == TRUE){
  d_habpref <- read.csv('~/secsse/habitat_preference_asv_new.csv')
  d_taxa <- readRDS('~/secsse/ps_otu_asv_filt.rds')
  tree <- read.tree('~/secsse/myxo_asv_chronopl10.tre')
  fit_musse_no_se <- readRDS('~/secsse/asv_musse_no_se.rds')
}

if(server == FALSE){
  # read in habitat preference
  d_habpref <- read.csv(here('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_asv_new.csv'))
  # read in phyloseq object and grab tax table
  d_taxa <- readRDS(here('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_asv_filt.rds'))
  # read in tree
  tree <- read.tree(here('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_chronopl10.tre'))
  # read in musse model
  fit_musse_no_se <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_musse_no_se.rds')
}

d_taxa <- d_taxa %>%
  phyloseq::tax_table() %>%
  data.frame() %>%
  janitor::clean_names() %>%
  rownames_to_column('otu')

# create d_meta
# use habitat_preference3 which collapses into marine mud generalists
d_meta <- left_join(select(d_habpref, otu, habitat_preference = habitat_preference3, num_present), select(d_taxa, otu:family))

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


# read in musse model
fit_musse_no_se <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_musse_no_se.rds')
#fit_musse_no_se <- readRDS('~/secsse/asv_musse_no_se.rds')

# try and run SecSSE which runs concealed state and speciation models
# https://cran.r-project.org/web/packages/secsse/vignettes/Using_secsse.html

# set up secsse model
trait <- data.frame(otu = names(hab_pref_num), trait = hab_pref_num)
traits <- sortingtraits(trait, tree)

# setup arguments to pass to secsse_ml

# set number of concealed states
num_concealed_states <- 2

# setup parameter list
idparslist <- id_paramPos(traits, num_concealed_states = num_concealed_states)

idparslist

# setup speciation rates ####
# first make all speciation rates the same within hidden states
idparslist$lambdas[] <- rep(1:2, each = 5)

# setup extinction rates ####
# firstly make all extinction rates the same
idparslist$mus[] <- 3

# setup transition rates ####

# make a bunch of transitions 0 
# these transitions were not possible in the markov model
# q14~0, q24~0, q25~0, q41~0, q42~0, q53~0, q45~0, q54~0
# with an ARD model cannot get away from hidden state models having more parameters!

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
         id = id*1000) 

# first make any of the transitions not possible in the Markov model 0
# q14~0, q24~0, q25~0, q41~0, q42~0, q53~0, q45~0, q54~0
q_matrix <- mutate(q_matrix,
                   new_id = ifelse(transition %in% c('q14', 'q24', 'q25', 'q41', 'q42', 'q53', 'q45', 'q54'), 0, id))

# make transitions that are across hidden states AND trait states 0. i.e. 1A -> 2B
q_matrix <- mutate(q_matrix,
                   new_id = ifelse(from_hidden != to_hidden & from_trait != to_trait, 0, new_id))

# make transitions that are across hidden states AND trait states 0. i.e. 1A -> 2B
q_matrix <- mutate(q_matrix,
                   new_id = ifelse(from_hidden != to_hidden & from_trait != to_trait, 0, new_id))


# work out which parameters need to be the same
# when the hidden transition is the same, give them the same value
hidden_to_assign <- select(q_matrix, from_hidden, to_hidden) %>%
  distinct() %>%
  mutate(hidden_id = letters[1:n()])

# when the trait transition is the same, give them the same value
trait_to_assign <- select(q_matrix, from_trait, to_trait) %>%
  distinct() %>%
  mutate(trait_id = 1:n())

# merge altogether with the q matrix
q_matrix <- left_join(q_matrix, hidden_to_assign) %>%
  left_join(trait_to_assign)

# if the trait state changes, is not 0 or NA, create a column giving it the trait id
q_matrix <- mutate(q_matrix, new_id2 = ifelse(from_trait != to_trait & new_id != 0 & !is.na(new_id), trait_id, new_id))
# if the hidden state changes, is not 0 or NA, create a column giving it the hidden id
q_matrix <- mutate(q_matrix, new_id3 = ifelse(from_hidden != to_hidden & new_id != 0 & !is.na(new_id), hidden_id, new_id))

# when numbers and letters are not in the possible set (from trait and hidden ID) (make them 0)
q_matrix <- mutate(q_matrix, new_id2 = ifelse(!new_id2 %in% trait_to_assign$trait_id, 0, new_id2),
                    new_id3 = ifelse(!new_id3 %in% hidden_to_assign$hidden_id, 0, new_id3))

# look at unique values of new id
unique(q_matrix$new_id3)

# make the final new ID
# make the numbers higher than would be possible for relabelling later
q_matrix <- mutate(q_matrix, new_id_final = case_when(new_id2 > 0 ~ as.character(new_id2),
                                                        new_id3 != 0 ~ new_id3,
                                                        TRUE ~ '0'),
                    new_id_final = as.numeric(as.factor(new_id_final)),
                    new_id_final = ifelse(new_id_final == 1, 0, new_id_final),
                    new_id_final = new_id_final*100)

# run a for loop to replace each number in the initial q matrix
q <- idparslist[[3]]

for(i in min(q,na.rm = TRUE):max(q, na.rm = TRUE)){
  q[which(q == i)] <- filter(q_matrix, (id/1000) == i) %>% pull(new_id_final)
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

# setup parameters that are not being estimated ####

# set q51 to be 3.5 times q15 to stop it going mad high
fit_musse_no_se$par
8.1/2.3

# set q15 to be 1.46 times q51 to stop it going mad high
fit_musse_no_se$par
11.87/8.12

# change q15 to be 3.5x q15
constraint_factors <- NULL # factors included in constraint
to_constrain <- c(15) # parameters to be a function of others q5a1a and q5b1b
init_factor <- NULL

# run secsse
functions_defining_params <- list()
functions_defining_params[[1]] <- function(){
  par_15 <- 1.46*par_6
}

# set initial values ####

# musse transition rates
musse_transitions <- fit_musse_no_se$par[grepl('q', names(fit_musse_no_se$par))]

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
  left_join(., data.frame(transition = names(musse_transitions), rate = unname(musse_transitions))) %>%
  mutate(rate = replace_na(rate, mean(rate, na.rm = TRUE))) %>%
  select(id, rate) %>%
  distinct() %>%
  filter(! id %in% to_constrain) %>%
  pull(rate)

initparsopt <- c(rep(init_lambda, times = 2),
                 rep(init_mu, times = 1),
                 init_transition)

# check number of estimated parameters is the same as number of initial values
idparsopt <- c(1:max(idparslist$Q, na.rm=TRUE))
idparsopt <- idparsopt[!idparsopt %in% to_constrain]

length(initparsopt) == length(idparsopt)

idparslist

# set the ID and values for the fixed parameters
idparsfix <- 0 # zeroes have the value of zero
parsfix <- 0

# set number of iterations
max_iter <- 1000 * round((1.25)^length(idparsopt))

# save out number length of start values as number of parameters
saveRDS(length(idparsopt), paste('~/secsse/nparams_', name, '.rds', sep =''))

# right think I have done it! Ridiculous
mod_secsse <- secsse_ml_func_def_pars(tree,
                                      traits, 
                                      num_concealed_states=num_concealed_states, 
                                      idparslist, 
                                      idparsopt, 
                                      initparsopt, 
                                      idparsfix, 
                                      parsfix, 
                                      cond="maddison_cond",
                                      root_state_weight = "maddison_weights", 
                                      tol = c(1e-04, 1e-05, 1e-07), 
                                      sampling_fraction=rep(1, times = length(unique(traits))), 
                                      maxiter = max_iter,
                                      use_fortran=TRUE,
                                      methode="ode45", 
                                      optimmethod = "simplex", 
                                      num_cycles = 1, 
                                      run_parallel=TRUE,
                                      idparsfuncdefpar = to_constrain,
                                      idfactorsopt = constraint_factors,
                                      initfactors = init_factor,
                                      functions_defining_params = functions_defining_params)

saveRDS(mod_secsse, paste('~/secsse/fit_', name, '.rds', sep =''))