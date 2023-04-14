# try a MuSSE using secSSE for the Sriswasdi data

# load in packages ####

# make sure curl is installed
library(curl)
librarian::shelf(diversitree, secsse, DDD, apTreeshape, doParallel, phytools, foreach, doMC, tidyverse, here, furrr)

# identify conflicts in the tidyverse packages and other packages
tidyverse_conflicts()

# load in data ####

# filename
name <- 'musse'

# server - yes or no
server <- FALSE

if(server == TRUE){
  # load in phylogenetic tree
  tree <- read.tree('sriswasdi/treecut_representative_topology_scale0.5_patched.newick')
  # load in trait data
  trait_file <- read.table('sriswasdi/treecut_representative_topology_scale0.5_patched_states.txt')
  # read in ARD Mk model
  fit_mk <- readRDS('sriswasdi/mod_ard.rds')
}

if(server == FALSE){
  # load in phylogenetic tree
  tree <- read.tree('data/sriswasdi_data/treecut_representative_topology_scale0.5_patched.newick')
  # load in trait data
  trait_file <- read.table('data/sriswasdi_data/treecut_representative_topology_scale0.5_patched_states.txt')
  # read in ARD Mk model
  fit_mk <- readRDS('data/sriswasdi_data/mod_ard.rds')
}

# add a very tiny number onto branch lengths that are zero
tree$edge.length[tree$edge.length == 0] <- 0.000001

# make tree ultrametric
tree <- phytools::force.ultrametric(tree, method = 'extend')

# turn trait into vector
trait <- trait_file[1,] |> as.numeric()
names(trait) <- colnames(trait_file)

# turn trait into dataframe
d_trait <- tibble(tip_label = names(trait),
                  trait_num = unname(trait),
                  trait = ifelse(trait_num == 0, 'specialist', 'generalist')) %>%
  left_join(tibble(tip_label = tree$tip.label), .)

head(tree$tip.label)

# check tip labels marry up to d_trait tip labels
sum(d_trait$tip_label == tree$tip.label) == length(tree$tip.label)
# SUCCESS if TRUE

# try and run SecSSE which runs concealed state and speciation models
# https://cran.r-project.org/web/packages/secsse/vignettes/Using_secsse.html

# set up secsse model
traits <- sortingtraits(select(d_trait, tip_label, trait_num), tree)

# setup arguments to pass to secsse_ml

# set number of concealed states
num_concealed_states <- 1

# setup parameter list
idparslist <- id_paramPos(traits, num_concealed_states = num_concealed_states)

idparslist

# setup speciation rates ####
# first make all speciation rates the same within hidden states
idparslist$lambdas[] <- 1:2

# setup extinction rates ####
# firstly make all extinction rates the same
idparslist$mus[] <- 3:4

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
         id = -id) # make id negative to prevent mismatches in ids when setting rules

colnames(q_matrix)

# make transitions that are across hidden states AND trait states 0. i.e. 1A -> 2B
q_matrix <- mutate(q_matrix,
                   new_id = ifelse(from_hidden != to_hidden & from_trait != to_trait, 0, id))

# work out which parameters need to be the same
# when the hidden transition is the same, give them the same value
hidden_to_assign <- select(q_matrix, from_hidden, to_hidden) %>%
  distinct() %>%
  mutate(hidden_id = letters[1:n()])

# give all trait transitions different values 1A -> 2A != 1B -> 2B
# do this by adding from_hidden to the select argument
trait_to_assign <- select(q_matrix, from_trait, to_trait, from_hidden) %>%
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
                                                      new_id3 != 0 ~ as.character(new_id3),
                                                      TRUE ~ '0'),
                   new_id_final = as.numeric(as.factor(new_id_final)),
                   new_id_final = ifelse(new_id_final == 1, 0, new_id_final),
                   new_id_final = new_id_final*100)

# run a for loop to replace each number in the initial q matrix
q <- idparslist[[3]]

for(i in min(q,na.rm = TRUE):max(q, na.rm = TRUE)){
  q[which(q == i)] <- filter(q_matrix, -id == i) %>% pull(new_id_final)
}

q

idparslist$Q <- q

# replace all values in q matrix that are not 0 or NA by next logical value

# find the value for the extinction parameter
mu <- max(idparslist$mus)

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
  distinct() %>%
  pull(rate)

initparsopt <- c(rep(init_lambda, times = max(idparslist$lambdas)),
                 rep(init_mu, times = length(unique(idparslist$lambdas))),
                 init_transition)

# check number of estimated parameters is the same as number of initial values
idparsopt <- c(1:max(idparslist$Q, na.rm=TRUE))

length(initparsopt) == length(idparsopt)

idparslist

# set the ID and values for the fixed parameters
idparsfix <- 0 # zeroes have the value of zero
parsfix <- 0

# set number of iterations
max_iter <- 1000 * round((1.25)^length(idparsopt))

# test initial values to work out their initial log likelihood
idparslist

# write function to get initial values into the correct format
get_inits_matrix <- function(inits, idparslist){
  for(i in 1:length(idparslist)){
    for(j in 1:length(idparslist[[i]])){
      if(idparslist[[i]][j] %in% 1:length(inits)){idparslist[[i]][j] <- inits[idparslist[[i]][j]]}
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
vals <- expand.grid(lambda = vals, mu = vals, q = vals) %>%
  as.tibble()

# create an empty dataframe
inits_ml <- mutate(vals, inits = list(NA), loglik = NA,
                   id = 1:n())

# set up for loop to run screen for initial values
pb <- progress::progress_bar$new(total = nrow(inits_ml))

for(i in 1:nrow(inits_ml)){
  pb$tick()
  
  # create inits
  temp_inits <- c((inits_lambda*inits_ml$lambda[i]),
                  (inits_mu*inits_ml$mu[i]),
                  (inits_q*inits_ml$q[i]))
  
  temp_mat <- get_inits_matrix(temp_inits, idparslist)
  
  temp_ml <- # check maximum likelihood values of initial values
    secsse_loglik(
      temp_mat,
      tree,
      traits,
      num_concealed_states = num_concealed_states,
      parsfix,
      cond = "maddison_cond",
      root_state_weight = "maddison_weights",
      sampling_fraction = sampled_fraction_1,
      see_ancestral_states = FALSE
    )
  
  inits_ml$inits[[i]] <- temp_inits
  
  inits_ml$loglik[i] <- temp_ml
  
}

# save this out
saveRDS(inits_ml, paste('data/sriswasdi_data/init_vals_ml/', name, '.rds', sep = ''))
