# script to run muhisse

# load in packages
library(hisse)
library(ape)

# load in necessary data
tree <- read.tree('data/sequencing_rpoB/processed/muhisse/tree_clean.tree')

# load in data for muhisse
muhisse_df <- readRDS('data/sequencing_rpoB/processed/muhisse/muhisse_trait_data.rds')

# load in coding to know what the 0s and 1s mean
coding <- readRDS('data/sequencing_rpoB/processed/muhisse/coding_for_states.rds')
# non-saline generalist is 00
# non-saline specialist is 01
# saline generalist is 10
# saline specialist is 11

# check the tip labels match the order of the otu names in the dataframe
sum(muhisse_df$tip_label == tree$tip.label) == length(tree$tip.label)
# TRUE is good

# create transition rate matrix
trans_mat <- TransMatMakerMuHiSSE(hidden.traits = 1)
trans_mat_nohidden <- TransMatMakerMuHiSSE(hidden.traits = 0)

# we now need to customise some of the transition rates and make some of them possible that are not currently
# namely 10 (saline generalism) to 01 (non-saline specialist)
# based on the Markov model that we did pre-musse
trans_mat['(10A)', '(01A)'] <- max(trans_mat, na.rm = TRUE) + 1
trans_mat['(10B)', '(01B)'] <- max(trans_mat, na.rm = TRUE) + 1

trans_mat_nohidden['(10)', '(01)'] <- max(trans_mat_nohidden, na.rm = TRUE) + 1

# As the best MuSSE model just has variation in speciation, this is what we will fit here
turnover <- c(1,2,3,4,5,6,7,8)
extinction_fraction <- rep(1, 8) 

turnover_nohidden <- c(1,2,3,4)
extinction_fraction_nohidden <- rep(1, 4)

# set estimated proportion of extant species - 1 for each
f = c(1,1,1,1)

# run muhisse with no hidden states
muhisse_nohidden <- MuHiSSE(phy=tree, data=muhisse_df, f=f,
                            turnover=turnover_nohidden, 
                            eps=extinction_fraction_nohidden, 
                            hidden.states=FALSE, 
                            trans.rate=trans_mat_nohidden)

saveRDS(muhisse_nohidden, 'data/sequencing_rpoB/processed/muhisse/muhisse_nohidden.rds')

# muhisse with 1 hidden state
muhisse_1 <- MuHiSSE(phy=tree, data=muhisse_df, f=f,
                     turnover=turnover, 
                     eps=extinction_fraction, 
                     hidden.states=TRUE, 
                     trans.rate=trans_mat)

saveRDS(muhisse_1, 'data/sequencing_rpoB/processed/muhisse/muhisse_1.rds')