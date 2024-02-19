# fit the MuSSE model using secsse

# load in packages ####

# make sure curl is installed
library(curl)
librarian::shelf(diversitree, rsetienne/secsse, DDD, doParallel, foreach, doMC, tidyverse, here, furrr)

# identify conflicts in the tidyverse packages and other packages
tidyverse_conflicts()

# function for running secsse
# write a custom function to do everything we want in terms of fitting the model and saving it out
fit_secsse <- function(list_inits_sampfrac){
  
  # pick out inits
  temp_inits <- list_inits_sampfrac$inits
  
  # pick out sampled_fractions
  temp_samp_frac <- list_inits_sampfrac$sampled_fractions
  
  # run secsse
  # right think I have done it! Ridiculous
  mod_secsse <- secsse_ml(
    tree,
    traits,
    num_concealed_states = num_concealed_states,
    idparslist,
    idparsopt,
    initparsopt = temp_inits,
    idparsfix,
    parsfix,
    cond = "maddison_cond",
    root_state_weight = "maddison_weights",
    tol = c(1e-04, 1e-05, 1e-07),
    sampling_fraction = temp_samp_frac,
    maxiter = max_iter,
    optimmethod = "simplex",
    num_cycles = 75,
    num_threads = 2,
    method = 'odeint::bulirsch_stoer',
    loglik_penalty = 0.1
  )
  
  # create a list of the output
  output <- list(n_params = length(idparsopt),
                 inits = temp_inits,
                 samp_frac = temp_samp_frac,
                 setup = idparslist,
                 mod = mod_secsse,
                 lambda_inits_fac = list_inits_sampfrac$lambda,
                 mu_inits_fac = list_inits_sampfrac$mu,
                 q_inits_fac = list_inits_sampfrac$q)
  
  # save out the list
  temp_name <- paste(name, '_', 'sampfrac', unique(temp_samp_frac), '_', 'run', list_inits_sampfrac$run,  sep = '')
  
  saveRDS(output, paste('~/secsse/seccse_', temp_name,  '_asvboot', num_tree, '.rds', sep =''))
  
}

# load in data ####

# filename
name <- 'musse'

# server - yes or no
server <- TRUE

# list bootstrapped trees to loop over
if(server == FALSE){
  trees <- list.files('data/sequencing_rpoB/raxml/trees/myxo_asv/bootstraps', full.names = TRUE, pattern = 'cv.tre$')
  mk_models <- list.files('data/sequencing_rpoB/processed/transition_rates', full.names = TRUE, pattern = 'asvboot')
  start_vals_files <- list.files('data/sequencing_rpoB/processed/secsse/init_vals_ml/bootstraps', full.names = TRUE)
}

if(server == TRUE){
  trees <- list.files('~/secsse/bootstraps/trees', full.names = TRUE, pattern = 'cv.tre$')
  mk_models <- list.files('~/secsse/bootstraps/mk_models', full.names = TRUE)
  start_vals_files <- list.files('~/secsse/bootstraps/start_vals', full.names = TRUE)
}


for(j in 1:length(trees)){
  
  if(server == FALSE){
    # read in habitat preference
    d_habpref <- read.csv(here('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_asv_new.csv'))
    # read in phyloseq object and grab tax table
    d_taxa <- readRDS(here('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_asv_filt.rds'))
    
    # read in tree
    tree <- read.tree(trees[j])
    
    num_tree <- parse_number(tools::file_path_sans_ext(basename(trees[j])))
    
    temp_mk <- mk_models[grepl(paste('asvboot', num_tree, sep =''), mk_models)]
    
    # read in Mk model
    fit_mk <- readRDS(temp_mk)
    
    start_vals <- readRDS(start_vals_files[grepl(paste(name, '_', 'asvboot', num_tree, sep =''), start_vals_files)])
  }
  
  
  if(server == TRUE){
    d_habpref <- read.csv('~/secsse/habitat_preference_asv_new.csv')
    d_taxa <- readRDS('~/secsse/ps_otu_asv_filt.rds')
    
    # read in tree
    tree <- read.tree(trees[j])
    
    num_tree <- parse_number(tools::file_path_sans_ext(basename(trees[j])))
    
    temp_mk <- mk_models[grepl(paste('asvboot', num_tree, sep =''), mk_models)]
    
    # read in Mk model
    fit_mk <- readRDS(temp_mk)
    
    # start vals
    start_vals <- readRDS(start_vals_files[grepl(paste(name, '_', 'asvboot', num_tree, sep =''), start_vals_files)])
    
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
  
  # setup different inits ####
  
  start_vals
  
  # filter out NaN and Inf
  start_vals <- filter(start_vals, !is.nan(loglik) & !is.infinite(loglik) & !is.na(loglik)) %>%
    # keep the six best log liks
    slice_max(., order_by = loglik, n = 6)
  
  inits <- start_vals$inits
  
  # also need to change the sample fractions ####
  # use 1, 0.5, 0.25, 0.125, 0.0625 to see how they change the fit
  
  sampled_fraction_1 <- rep(1, times = length(unique(traits)))
  sampled_fraction_0.5 <- rep(0.5, times = length(unique(traits)))
  sampled_fraction_0.25 <- rep(0.25, times = length(unique(traits)))
  sampled_fraction_0.125 <- rep(0.125, times = length(unique(traits)))
  sampled_fraction_0.0625 <- rep(0.0625, times = length(unique(traits)))
  
  sampled_fractions <- list(sampled_fraction_0.0625, sampled_fraction_0.125, sampled_fraction_0.25, sampled_fraction_0.5, sampled_fraction_1)
  
  num_samp_frac <- length(sampled_fractions)
  
  # create all combinations of the two lists
  all_combs <- expand_grid(sampled_fractions, inits) %>%
    mutate(run = 1:n(),
           lambda = rep(start_vals$lambda, times = num_samp_frac),
           mu = rep(start_vals$mu, times = num_samp_frac)) %>%
    purrr::transpose()
  
  # Set a "plan" for how the code should run.
  plan(multisession, workers = 10)
  
  # run future_walk
  furrr::future_walk(all_combs, fit_secsse)
  
}

# set a test
secsse_ml(
  tree,
  traits,
  num_concealed_states = num_concealed_states,
  idparslist,
  idparsopt,
  initparsopt = all_combs[[1]]$inits,
  idparsfix,
  parsfix,
  cond = "maddison_cond",
  root_state_weight = "maddison_weights",
  tol = c(1e-04, 1e-05, 1e-07),
  sampling_fraction = all_combs[[1]]$sampled_fractions,
  maxiter = max_iter,
  optimmethod = "simplex",
  num_cycles = 75,
  num_threads = 2,
  method = 'odeint::runge_kutta_cash_karp54',
  loglik_penalty = 0.05
)
