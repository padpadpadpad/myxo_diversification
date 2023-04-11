# try and recreate the bisse analysis of Sriswasdi et al 2016

# load in packages
librarian::shelf(diversitree, ape, ggtree, phytools, tidyverse, hisse)

tidyverse::tidyverse_conflicts()

# load in phylogenetic tree
tree <- read.tree('data/sriswasdi_data/treecut_representative_topology_scale0.5_patched.newick')
tree

is.ultrametric(tree)

tree <- phytools::force.ultrametric(tree, method = 'extend')

ggtree(tree, layout = 'circular')

# read in trait data
# 399 generalists and 223
# from the paper
# so 0 is specialist, 1 is generalist
trait_file <- read.table('data/sriswasdi_data/treecut_representative_topology_scale0.5_patched_states.txt')

str(trait_file)

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

# plot tree and traits as colours
ggtree(tree, layout = 'dendrogram') %<+% d_trait +
  geom_tippoint(aes(x=x+0.1, col = trait), size = 0.4) +
  guides(color = guide_legend(override.aes = list(size = 5)))

# ok fit a Mk model with all rates different and symmetrical/equal rates using diversitree

# set up models
lik_ard <- make.mk2(tree, trait)
lik_sym <- constrain(lik_ard, q01~q10)

argnames(lik_ard)
argnames(lik_sym)

# set initial values
inits_ard <- rep(1, length(argnames(lik_ard)))
inits_sym <- rep(1, length(argnames(lik_sym)))

mod_ard <- find.mle(lik_ard, inits_ard, method = 'subplex', control = list(maxit = 50000))
mod_sym <- find.mle(lik_sym, inits_sym, method = 'subplex', control = list(maxit = 50000))

AIC(mod_ard, mod_sym) %>% data.frame() %>%
  mutate(weights = MuMIn::Weights(AIC))
# mod ard is favoured

mod_ard$par

# setup bisse ####

# set up sampling fractions, set them all to 1
sampling_frac <- setNames(rep(1, times = 2), sort(unique(trait)))

# set up likelihood model for diversitree
lik_bisse <- make.bisse(tree, trait, sampling.f = sampling_frac)

# we can estimate starting values using starting.point.bisse()
start_vals <- starting.point.bisse(tree)

# replace the start values with those in the best Mk model
for(i in 1:length(mod_ard$par)){
  start_vals[names(start_vals) == names(mod_ard$par)[i]] <- unname(mod_ard$par[i])
}

# fit musse model
fit_bisse <- find.mle(lik_bisse, x.init = start_vals[argnames(lik_bisse)], method = 'subplex', control = list(maxit = 50000))
fit_bisse$par
fit_bisse1 <- fit_bisse

# remove state dependent speciation and extinction parameters
lik_null_no_sse <- constrain(lik_bisse, 
                             lambda1 ~ lambda0,
                             mu1 ~ mu0)

# remove only speciation parameters
lik_null_no_ss <- constrain(lik_bisse, 
                            lambda1 ~ lambda0)

# remove only extinction parameters
lik_null_no_se <- constrain(lik_bisse, mu1 ~ mu0)

# fit model
fit_bisse_no_sse <- find.mle(lik_null_no_sse, x.init = start_vals[argnames(lik_null_no_sse)], method = 'subplex', control = list(maxit = 100000))
fit_bisse_no_ss <- find.mle(lik_null_no_ss, x.init = start_vals[argnames(lik_null_no_ss)],method = 'subplex', control = list(maxit = 50000))
fit_bisse_no_se <- find.mle(lik_null_no_se, x.init = start_vals[argnames(lik_null_no_se)], method = 'subplex', control = list(maxit = 100000))

# make AIC table
d_bisse_aic <- data.frame(model = c('full_sse', 'no_se', 'no_ss', 'no_sse'), 
                          aic = c(AIC(fit_bisse), AIC(fit_bisse_no_se), AIC(fit_bisse_no_ss), AIC(fit_bisse_no_sse))) %>%
  dplyr::arrange(., aic) %>%
  mutate(weights = round(MuMIn::Weights(aic), 3))

d_bisse_aic

d_bisse_aic1 <- d_bisse_aic
# no se model is favoured - makes fitting the hisse easier

# run models again but with ultrametric tree with no zero branches
tree2 <- read.tree('data/sriswasdi_data/treecut_representative_topology_scale0.5_patched.newick')

# add a very tiny number onto branch lengths that are zero
tree2$edge.length[tree$edge.length == 0] <- 0.000001

# make tree ultrametric
tree2 <- phytools::force.ultrametric(tree2, method = 'extend')

# set up sampling fractions, set them all to 1
sampling_frac <- setNames(rep(1, times = 2), sort(unique(trait)))

# set up likelihood model for diversitree
lik_bisse <- make.bisse(tree2, trait, sampling.f = sampling_frac)

# we can estimate starting values using starting.point.bisse()
start_vals <- starting.point.bisse(tree)

# replace the start values with those in the best Mk model
for(i in 1:length(mod_ard$par)){
  start_vals[names(start_vals) == names(mod_ard$par)[i]] <- unname(mod_ard$par[i])
}

# fit musse model
fit_bisse <- find.mle(lik_bisse, x.init = start_vals[argnames(lik_bisse)], method = 'subplex', control = list(maxit = 50000))
fit_bisse$par

# remove state dependent speciation and extinction parameters
lik_null_no_sse <- constrain(lik_bisse, 
                             lambda1 ~ lambda0,
                             mu1 ~ mu0)

# remove only speciation parameters
lik_null_no_ss <- constrain(lik_bisse, 
                            lambda1 ~ lambda0)

# remove only extinction parameters
lik_null_no_se <- constrain(lik_bisse, mu1 ~ mu0)

# fit model
fit_bisse_no_sse <- find.mle(lik_null_no_sse, x.init = start_vals[argnames(lik_null_no_sse)], method = 'subplex', control = list(maxit = 100000))
fit_bisse_no_ss <- find.mle(lik_null_no_ss, x.init = start_vals[argnames(lik_null_no_ss)],method = 'subplex', control = list(maxit = 50000))
fit_bisse_no_se <- find.mle(lik_null_no_se, x.init = start_vals[argnames(lik_null_no_se)], method = 'subplex', control = list(maxit = 100000))

# make AIC table
d_bisse_aic2 <- data.frame(model = c('full_sse', 'no_se', 'no_ss', 'no_sse'), 
                          aic = c(AIC(fit_bisse), AIC(fit_bisse_no_se), AIC(fit_bisse_no_ss), AIC(fit_bisse_no_sse))) %>%
  dplyr::arrange(., aic) %>%
  mutate(weights = round(MuMIn::Weights(aic), 3))

d_bisse_aic1
d_bisse_aic2

# do not run this code ####

# set up the hisse model

# make a null model
null_rate_matrix <- TransMatMakerHiSSE(hidden.traits = 1, make.null=TRUE)
null_rate_matrix

# The states are ordered 0A, 1A, 0B, 1B, so we will specify that the diversification rates are identical between 0 and 1 but can differ between A and B - i.e. that states 0A and 1A share one value, and 0B and 1B share another. 
# i.e. changes in diversification are solely due to the hidden state
null_net_turnover <- c(1,1,2,2)
null_extinction_fraction <- c(1,1,2,2)

# make dataframe for passing to hisse
d_hisse <- select(d_trait, tip_label, trait_num) %>% data.frame()
row.names(d_hisse) <- d_hisse$tip_label

# run null hisse
null_hisse <- hisse(phy = tree, 
                    data = d_hisse, 
                    hidden.states=TRUE,
                    turnover = null_net_turnover, 
                    eps = null_extinction_fraction,
                    trans.rate = null_rate_matrix)

saveRDS(mod_ard, 'data/sriswasdi_data/mod_ard.rds')
