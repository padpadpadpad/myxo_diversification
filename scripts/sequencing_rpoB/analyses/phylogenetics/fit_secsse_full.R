# try a multistate MuHiSSE using secSSE


# load in packages ####

librarian::shelf(here, ggtree, diversitree, tidybayes, ggdist, secsse, DDD, tidyverse, apTreeshape, doParallel, foreach, doMC)

# identify conflicts in the tidyverse packages and other packages
tidyverse_conflicts()

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

# read in musse model
fit_musse_no_se <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_musse_no_se.rds')

# try and run SecSSE which runs concealed state and speciation models
# https://cran.r-project.org/web/packages/secsse/vignettes/Using_secsse.html

trait <- data.frame(otu = names(hab_pref_num), trait = hab_pref_num)
traits <- sortingtraits(trait, tree)

# setup arguments to pass to secsse_ml

num_concealed_states <- 2

# setup IDPARSLIST
idparslist <- id_paramPos(traits, num_concealed_states = num_concealed_states)

idparslist

# firstly make all extinction rates the same
idparslist$mus[] <- 11

# make a bunch of transitions 0 

# these transitions were not possible in the markov model
# q14~0, q24~0, q25~0, q41~0, q42~0, q53~0, q45~0, q54~0

# first set all of these transitions to zero, including the hidden state versions
idparslist[[3]][1, c(4, 9)] <- 0 # q14
idparslist[[3]][6, c(4, 9)] <- 0 # q14
idparslist[[3]][2, c(4, 9)] <- 0 # q24
idparslist[[3]][7, c(4, 9)] <- 0 # q24
idparslist[[3]][5, c(4, 9)] <- 0 # q54
idparslist[[3]][10, c(4, 9)] <- 0 # q54
idparslist[[3]][4, c(1, 6)] <- 0 # q41
idparslist[[3]][9, c(1, 6)] <- 0 # q41
idparslist[[3]][4, c(2, 7)] <- 0 # q42
idparslist[[3]][9, c(2, 7)] <- 0 # q42
idparslist[[3]][4, c(5, 10)] <- 0 # q45
idparslist[[3]][9, c(5, 10)] <- 0 # q45
idparslist[[3]][2, c(5, 10)] <- 0 # q25
idparslist[[3]][7, c(5, 10)] <- 0 # q25
idparslist[[3]][5, c(3, 8)] <- 0 # q53
idparslist[[3]][10, c(3, 8)] <- 0 # q53

idparslist

# now disallow dual transitions - cannot go from 1A to 2B etc
idparslist[[3]][1, c(7,8,9,10)] <- 0
idparslist[[3]][2, c(6,8,9,10)] <- 0
idparslist[[3]][3, c(6,7,9,10)] <- 0
idparslist[[3]][4, c(6,7,8,10)] <- 0
idparslist[[3]][5, c(6,7,8,9)] <- 0
idparslist[[3]][6, c(2,3,4,5)] <- 0
idparslist[[3]][7, c(1,3,4,5)] <- 0
idparslist[[3]][8, c(1,2,4,5)] <- 0
idparslist[[3]][9, c(1,2,3,5)] <- 0
idparslist[[3]][10, c(1,2,3,4)] <- 0

idparslist

# replace all values in q matrix that are not 0 or NA by next logical value
num_q_params <- length(idparslist$Q[!is.na(idparslist$Q) & idparslist$Q != 0])
idparslist$Q[!is.na(idparslist$Q) & idparslist$Q != 0] <- 12:(12+num_q_params-1)

num_params_musse <- length(fit_musse_no_se$par)

# set q51 to be 3.5 times q15 to stop it going mad high
fit_musse_no_se$par
8.1/2.3

# change q5b1b to be 3.5x q1b5b
# change q5a1a to be 3.5x q1a5a
constraint_factors <- NULL # factors included in constraint
to_constrain <- c(14, 32) # parameters to be a function of others q5a1a and q5b1b
init_factor <- NULL

# run secsse
functions_defining_params <- list()
functions_defining_params[[1]] <- function(){
  par_14 <- 3.5*par_26
}
functions_defining_params[[2]] <- function(){
  par_32 <- 3.5*par_44
}

# set initial values

# set initial values for transition rates based on those from diversitree
# remember no q51
first_hidden_state <- unname(c(fit_musse_no_se$par[7:9], 0.25, fit_musse_no_se$par[10:11], 0.25, fit_musse_no_se$par[12:15], 0.25, fit_musse_no_se$par[16], 0.25, fit_musse_no_se$par[18], 0.25))
second_hidden_state <- unname(c(0.25, fit_musse_no_se$par[7:9], 0.25, fit_musse_no_se$par[10:11], 0.25, fit_musse_no_se$par[12:15], 0.25, fit_musse_no_se$par[16], 0.25, fit_musse_no_se$par[18]))
initparsopt <- c(rep(unname(fit_musse_no_se$par[1:5]), times = 2),
                 rep(unname(fit_musse_no_se$par[6]), times = 1),
                 c(first_hidden_state, second_hidden_state))

idparsopt <- c(1:45)
idparsopt <- idparsopt[!idparsopt %in% to_constrain]

length(initparsopt) == length(idparsopt)

idparsfix <- 0
parsfix <- 0

max_iter <- 1000 * round((1.25)^length(idparsopt))

# right think I have done it! Ridiculous
mod_secsse <- secsse_ml_func_def_pars(tree,
          traits, 
          num_concealed_states=2, 
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



mod_secsse <- secsse_ml(tree,
          traits, 
          num_concealed_states=2, 
          idparslist, 
          idparsopt, 
          initparsopt, 
          idparsfix, 
          parsfix, 
          cond="maddison_cond",
          root_state_weight = "maddison_weights", 
          tol = c(1e-04, 1e-05, 1e-07), 
          sampling_fraction=rep(1, times = length(unique(traits))), 
          maxiter = 1000 * round((1.25)^length(idparsopt)), 
          use_fortran=TRUE,
          methode="ode45", 
          optimmethod = "simplex", 
          num_cycles = 1, 
          run_parallel=TRUE)
