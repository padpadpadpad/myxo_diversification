# check bamm run convergence

# load packages
library(BAMMtools)
library(coda)

# set file
file <- 'data/sequencing_rpoB/bamm/bamm_asv_mcmc_out.txt'

# read in mcmc output
mcmcout <- read.csv(file, header=TRUE)
max(mcmcout$generation)

# discard some runs as burnin. We will discard the first 10% of samples
burnstart <- floor(0.4 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

# calculate effective sample size
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
