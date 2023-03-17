# lets look at which hidden model is best

# load packages
librarian::shelf(secsse, tidyverse)

# load in models
model_two_trait <- readRDS('data/sequencing_rpoB/processed/secsse/fit_twohidden_traitspeciation.rds')
model_two_hidden <- readRDS('data/sequencing_rpoB/processed/secsse/fit_twohidden_notraitspeciation.rds')
model_two_traithidden <- readRDS('data/sequencing_rpoB/processed/secsse/fit_twohidden_traithiddenspeciation.rds')

model_two_trait$MLpars
model_two_hidden$MLpars
model_two_traithidden$MLpars

nparams_two_trait <- readRDS('data/sequencing_rpoB/processed/secsse/nparams_twohidden_traitspeciation.rds')
nparams_two_hidden <- readRDS('data/sequencing_rpoB/processed/secsse/nparams_twohidden_notraitspeciation.rds')
nparams_two_traithidden <- readRDS('data/sequencing_rpoB/processed/secsse/nparams_twohidden_traithiddenspeciation.rds')

# read in musse model
fit_musse_no_se <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_musse_no_se.rds')
AIC(fit_musse_no_se)

d_results <- tibble(
  nhidden = c(2,2,2,0),
  model = c('trait_only', 'hidden_only', 'trait_and_hidden', 'musse'),
  loglik = c(model_two_trait$ML, model_two_hidden$ML, model_two_traithidden$ML, fit_musse_no_se$lnLik),
  nparams = c(nparams_two_trait, nparams_two_hidden, nparams_two_traithidden, length(fit_musse_no_se$par)),
  aic = (2*nparams) - 2*loglik
)
