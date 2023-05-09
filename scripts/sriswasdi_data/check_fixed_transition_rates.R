# compare Sriswasdi musse with ARD and ER to try and limit the number of estimated parameters

# load in packages
librarian::shelf(diversitree, ape, ggtree, phytools, tidyverse, hisse)

tidyverse::tidyverse_conflicts()

# run models again but with ultrametric tree with no zero branches
tree <- read.tree('data/sriswasdi_data/treecut_representative_topology_scale0.5_patched.newick')

# add a very tiny number onto branch lengths that are zero
tree$edge.length[tree$edge.length == 0] <- 0.000001

# make tree ultrametric
tree <- phytools::force.ultrametric(tree, method = 'extend')

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

# fit model where we constrain values

# 1. make transition rates the same
lik_null_no_se2 <- constrain(lik_bisse, mu1 ~ mu0, q01 ~ q10)
fit_bisse_no_se2 <- find.mle(lik_null_no_se2, x.init = start_vals[argnames(lik_null_no_se2)], method = 'subplex', control = list(maxit = 100000))
lik_null_no_se3 <- constrain(lik_bisse, mu1 ~ mu0, 
                             q01 ~ 2.553931,
                             q10 ~ 4.938238)
fit_bisse_no_se3 <- find.mle(lik_null_no_se3, x.init = start_vals[argnames(lik_null_no_se3)], method = 'subplex', control = list(maxit = 100000))

data.frame(aic = c(AIC(fit_bisse_no_se), AIC(fit_bisse_no_se2), AIC(fit_bisse_no_se3))) %>%
  dplyr::arrange(., aic) %>%
  mutate(weights = round(MuMIn::Weights(aic), 3))

fit_bisse_no_se3$par
fit_bisse_no_se$par
