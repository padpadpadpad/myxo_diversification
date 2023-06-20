#----------------------------------#
# lets look at phylogenetic signal #
#----------------------------------#

# resources:
# https://nhcooper123.github.io/pcm-primer-online/phylogenetic-signal-in-r.html
# https://twitter.com/padpadpadpad/status/1484458056419315712
# https://github.com/mrborges23/delta_statistic

#--------------------------#
# what this script does ####
#--------------------------#

# read in every non-ultrametric, rooted tree
# estimates phylogenetic signal a few different ways
# estimates the D statistic but variable needs to be binary
# estimates Pagels lambda using fitDicrete()

# load packages
library(tidyverse)
library(ape)
library(geiger)
library(caper)
library(ggtree)
library(ape)
library(here)

# set variable for ASV or otu similarity
otu_similarity <- c(99:91, 97.7, 'asv')

# set number of iterations or permutations
n_iter = 50 # number of iterations for calculating Pagel's lambda
n_permute = 1000 # number of permutations for calculating the D statistic

# set number of boot straps for subsampling the tree
n_boots = 100

# set number of cores
n_cores = 4
options(mc.cores=4)

# set number of tips to subsample
# the 91% tree is 353 tips
n_tips = 350

# run for loop for each otu similarity
for(i in 1:length(otu_similarity)){
  
  # set otu similarity
  temp_otu_similarity = otu_similarity[i]
  
  # read in dataset about habitat preferences
  if(temp_otu_similarity == 'asv'){d_pref <- read.csv(here(paste('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_', temp_otu_similarity, '.csv', sep = ''))) %>%
    dplyr::select(otu, habitat_preference)}
  
  if(temp_otu_similarity != 'asv'){d_pref <- read.csv(here(paste('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_', temp_otu_similarity, 'percent.csv', sep = ''))) %>%
    dplyr::select(otu, habitat_preference)}
  
  # read in rooted tree - tree has been rooted in figtree
  tree <- ape::read.tree(paste('data/sequencing_rpoB/raxml/trees/myxo_', temp_otu_similarity,'/myxo_', temp_otu_similarity, '_chronopl10.tre', sep = ''))
  
  # reorder d_pref otu to be the same order as tip labels
  d_pref2 <- tibble(otu = tree$tip.label) %>%
    left_join(., d_pref)
  
  # create vector of habitat preference
  hab_pref_vec <- d_pref2$habitat_preference
  names(hab_pref_vec) <- d_pref2$otu
  
  #----------------------------------------------------------------#
  # method 1. lets try get estimate D to fit on a binary trait. ####
  #----------------------------------------------------------------#
  
  # can only calculate D on binary traits, so will make three new vectors
  # can you live in freshwater (or not), live in marine mud (or not), or live on land (or not)
  d_pref2 <- mutate(d_pref2, freshwater = ifelse(str_detect(habitat_preference, 'freshwater'), 1, 0),
                    marine_mud = ifelse(str_detect(habitat_preference, 'marine_mud'), 1, 0),
                    terrestrial = ifelse(str_detect(habitat_preference, 'terrestrial'), 1, 0),
                    across(freshwater:terrestrial, ~replace_na(.x, 0))) %>%
    data.frame()
  
  # create combined dataset of phylogenetic tree and habitat preference data
  comb_data <- comparative.data(phy = tree, data = d_pref2, 
                                names.col = 'otu', vcv = TRUE,
                                na.omit = FALSE, warn.dropped = TRUE)
  
  # check dropped species
  comb_data$dropped$unmatched.rows
  comb_data$dropped$tips
  
  # delete any 0 length branches and collapse them into polytomies
  comb_data$phy <- di2multi(comb_data$phy)
  
  # estimate d for freshwater
  d_freshwater <- phylo.d(data = comb_data, names.col = otu, binvar = freshwater, 
                          permut = n_permute)
  d_freshwater
  plot(d_freshwater)
  
  # estimate d for marine mud
  d_marine_mud <- phylo.d(data = comb_data, names.col = otu, binvar = marine_mud, 
                         permut = n_permute)
  d_marine_mud
  plot(d_marine_mud)
  
  # estimate d for the terrestrial environment
  d_terrestrial <- phylo.d(data = comb_data, names.col = otu, binvar = terrestrial,
                           permut = n_permute)
  d_terrestrial
  plot(d_terrestrial)
  
  #---------------------------------------------------------------------------#
  # method 2. fit model for discrete character evolution using fitDiscrete ####
  #---------------------------------------------------------------------------#
  
  # convert habitat preference vector into all the binary instances
  freshwater_vec <- d_pref2$freshwater
  freshwater_vec <- ifelse(freshwater_vec == 1, 'freshwater', 'non_freshwater')
  marine_mud_vec <- d_pref2$marine_mud
  marine_mud_vec <- ifelse(marine_mud_vec == 1, 'marine_mud', 'non_marine_mud')
  terrestrial_vec <- d_pref2$terrestrial
  terrestrial_vec <- ifelse(terrestrial_vec == 1, 'terrestrial', 'non_terrestrial')
  
  names(freshwater_vec) <- names(terrestrial_vec) <- names(marine_mud_vec) <- d_pref2$otu
  
  # first fit for all traits - use ER as SYM and ARD takes SO LONG
  # run it for less iterations
  mod_all <- fitDiscrete(phy = tree, dat = hab_pref_vec, model = 'ER', transform = 'lambda', ncores = n_cores, control = list(niter = 100))
  
  # second for all freshwater
  mod_freshwater <- fitDiscrete(phy = tree, dat = freshwater_vec, model = 'ARD', transform = 'lambda', ncores = n_cores, control = list(niter = n_iter))
  
  # mud and shore
  mod_mud <- fitDiscrete(phy = tree, dat = marine_mud_vec, model = 'ARD', transform = 'lambda', ncores = n_cores, control = list(niter = n_iter))
  
  # terrestrial
  mod_terrestrial <- fitDiscrete(phy = tree, dat = terrestrial_vec, model = 'ARD', transform = 'lambda', ncores = n_cores, control = list(niter = n_iter))
  
  #---------------------------#
  # make output dataframes ####
  #---------------------------#
  
  temp_d_statistics <- tibble(habitat_preference = c('freshwater', 'marine_mud', 'terrestrial'),
                              estimated_d = c(d_freshwater$DEstimate, d_marine_mud$DEstimate, d_terrestrial$DEstimate),
                              prob_random = c(d_freshwater$Pval1, d_marine_mud$Pval1, d_terrestrial$Pval1),
                              prob_brownian = c(d_freshwater$Pval0, d_marine_mud$Pval0, d_terrestrial$Pval0),
                              similarity = temp_otu_similarity)
  
  temp_lambda <- tibble(habitat_preference = rep(c('all', 'freshwater', 'marine_mud', 'terrestrial'), each = 1),
                        estimated_lambda = c(mod_all$opt$lambda, mod_freshwater$opt$lambda, mod_mud$opt$lambda, mod_terrestrial$opt$lambda),
                        model = c('er', 'ard', 'ard', 'ard'),
                        similarity = temp_otu_similarity)
  
  # save out output
  write.csv(temp_d_statistics, here(paste('data/sequencing_rpoB/processed/phylogenetic_signal/d_statistic_', temp_otu_similarity, '.csv', sep = '')))
  write.csv(temp_lambda, here(paste('data/sequencing_rpoB/processed/phylogenetic_signal/lambda_', temp_otu_similarity, '.csv', sep = '')))
  
  #---------------------------------------------------------------#
  # subsample tree to 1000 tips and calculate the same metrics ####
  #---------------------------------------------------------------#
  
  # set up empty dataframes to populate bootstrapped results
  
  # setup empty d statistic dataframe
  temp_boot_d_statistics <- tibble(boot_num = 1:n_boots) %>%
    group_by(boot_num) %>%
    tidyr::expand(habitat_preference = temp_d_statistics$habitat_preference) %>%
    mutate(estimated_d = NA,
           prob_random = NA,
           prob_brownian = NA) %>%
    ungroup() %>%
    mutate(similarity = temp_otu_similarity)
  
  # set up empty lambda dataframe
  temp_boot_lambda <- tibble(boot_num = 1:n_boots) %>%
    group_by(boot_num) %>%
    tidyr::expand(habitat_preference = rep(c('all', 'freshwater', 'marine_mud', 'terrestrial'), each = 1)) %>%
    mutate(model = c('er', 'ard', 'ard', 'ard')) %>%
    ungroup() %>%
    mutate(similarity = temp_otu_similarity,
           estimated_lambda = NA) 
  
  # run for loop
  for(j in 1:n_boots){
    
    # calculate rows which have boot in
    rows_d_stat <- which(temp_boot_d_statistics[,1] == j)
    rows_lambda <- which(temp_boot_lambda[,1] == j)
    
    # sample tips to keep
    keep_tips <- sample(tree$tip.label, n_tips)
    tree_sub <- keep.tip(tree, keep_tips)
    
    # calculate D statistic
    d_pref2_sub <- filter(d_pref2, otu %in% tree_sub$tip.label) %>%
      data.frame()
    
    # create combined dataset of phylogenetic tree and preference for mud and shore
    comb_data <- comparative.data(phy = tree_sub, data = d_pref2_sub, 
                                  names.col = 'otu', vcv = TRUE,
                                  na.omit = FALSE, warn.dropped = TRUE)
    
    # delete any 0 length branches and collapse them into polytomies
    comb_data$phy <- di2multi(comb_data$phy)
    
    # estimate d
    d_freshwater <- phylo.d(data = comb_data, names.col = otu, binvar = freshwater,
                            permut = n_permute)
    
    d_marine_mud <- phylo.d(data = comb_data, names.col = otu, binvar = marine_mud,
                           permut = n_permute)
    
    d_terrestrial <- phylo.d(data = comb_data, names.col = otu, binvar = terrestrial, 
                             permut = n_permute)
    
    # calculate lamba
    hab_pref_vec_sub <- hab_pref_vec[names(hab_pref_vec) %in% keep_tips]
    freshwater_vec_sub <- d_pref2_sub$freshwater
    freshwater_vec_sub <- ifelse(freshwater_vec_sub == 1, 'freshwater', 'non_freshwater')
    marine_mud_vec_sub <- d_pref2_sub$marine_mud
    marine_mud_vec_sub <- ifelse(marine_mud_vec_sub == 1, 'marine_mud', 'non_marine_mud')
    terrestrial_vec_sub <- d_pref2_sub$terrestrial
    terrestrial_vec_sub <- ifelse(terrestrial_vec_sub == 1, 'terrestrial', 'non_terrestrial')
    
    names(freshwater_vec_sub) <- names(terrestrial_vec_sub) <- names(marine_mud_vec_sub) <- d_pref2_sub$otu
    
    # first fit for all traits - Use ER because otherwise this step is going take forever - especially for the bootstraps
    mod_all <- fitDiscrete(phy = tree_sub, dat = hab_pref_vec_sub, model = 'ER', transform = 'lambda', ncores = n_cores, control = list(niter = 20))
    
    # second for all freshwater
    mod_freshwater <- fitDiscrete(phy = tree_sub, dat = freshwater_vec_sub, model = 'ARD', transform = 'lambda', ncores = n_cores, control = list(niter = n_iter))
    
    # mud and shore
    mod_mud <- fitDiscrete(phy = tree_sub, dat = marine_mud_vec_sub, model = 'ARD', transform = 'lambda', ncores = n_cores, control = list(niter = n_iter))
    
    # terrestrial
    mod_terrestrial <- fitDiscrete(phy = tree_sub, dat = terrestrial_vec_sub, model = 'ARD', transform = 'lambda', ncores = n_cores, control = list(niter = n_iter))
    
    # attach info to precreated dataset
    temp_boot_d_statistics$estimated_d[rows_d_stat] <-c(d_freshwater$DEstimate, d_marine_mud$DEstimate, d_terrestrial$DEstimate)
    temp_boot_d_statistics$prob_random[rows_d_stat] <- c(d_freshwater$Pval1, d_marine_mud$Pval1, d_terrestrial$Pval1)
    temp_boot_d_statistics$prob_brownian[rows_d_stat] <- c(d_freshwater$Pval0, d_marine_mud$Pval0, d_terrestrial$Pval0)
    
    temp_boot_lambda$estimated_lambda[rows_lambda] <- c(mod_all$opt$lambda, mod_freshwater$opt$lambda,  mod_mud$opt$lambda, mod_terrestrial$opt$lambda)
  }

# save out dataset
write.csv(temp_boot_d_statistics, here(paste('data/sequencing_rpoB/processed/phylogenetic_signal/d_statistic_boot_', temp_otu_similarity, '.csv', sep = '')))
write.csv(temp_boot_lambda, here(paste('data/sequencing_rpoB/processed/phylogenetic_signal/lambda_boot_', temp_otu_similarity, '.csv', sep = '')))

# remove big objects
rm(list = c('d_terrestrial', 'd_marine_mud', 'd_freshwater', 'comb_data', 'mod_all', 'mod_mud', 'mod_freshwater', 'mod_terrestrial', 'tree_sub', 'd_pref2_sub', 'temp_boot_d_statistics', 'temp_boot_lambda', 'terrestrial_vec_sub', 'hab_pref_vec_sub', 'freshwater_vec_sub', 'marine_mud_vec_sub'))
  
}
