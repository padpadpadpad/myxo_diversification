# lets see if we can play around with the stochastic character mapping

library(phytools)
library(progress)
library(here)
library(tidyverse)
library(castor)
library(tidytree)


# reload simmap files in
simmap_files <- list.files(here("data/sequencing_rpoB/processed/transition_rates/simmap"), full.names = TRUE)
simmap_files <- simmap_files[simmap_files != here('data/sequencing_rpoB/processed/transition_rates/simmap/simmap_summary.rds')]

simmap_best <- purrr::map(simmap_files, readRDS)

simmap_best <- do.call(c, simmap_best)

# load in summary as well!
simmap_summary <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/simmap/simmap_summary.rds'))
# try and write an example to grab the mapped edges of all transitions
# grab out code from countSimmap()

tree <- simmap_best[[1]]

bb <- function(tree){
  # calculate total number of transitions
  n <- sum(sapply(tree$maps, length)) - nrow(tree$edge)
  
  # set states
  states <- colnames(tree$mapped.edge)
  
  # length of states
  m <- length(states)
  
  TT <- matrix(NA, m, m, dimnames = list(states, states))
  
  # write function to count transitions
  gg <- function(map, a, b) {
    
    if (length(map) == 1) 
      zz <- 0
    else {
      zz <- 0
      i <- 2
      while (i <= length(map)) {
        if (names(map)[i] == b && names(map)[i - 1] == 
            a) 
          zz <- zz + 1
        i <- i + 1
      }
    }
    return(zz)
  }
  
  
  for (i in 1:m) for (j in 1:m){
    if (i == j) TT[i, j] <- 0
    else TT[i, j] <- sum(sapply(tree$maps, gg, a = states[i], 
                                b = states[j]))
  }
  
  TT
  
  ###
  # New code to extract edge position of each transition ####
  ###
  
  # find positions in matrix that arent 0
  non_zero <- which(TT>0, arr.ind=TRUE)
  non_zero <- data.frame(non_zero)
  
  # create an empty list of this length
  res <- vector(mode='list', length=nrow(non_zero))
  
  # write a function to return mapped edges where transition occurred
  aa <- function(tree, a, b){
    temp_map <- sapply(tree$maps, gg, a = a, 
                       b = b)
    
    to_return <- data.frame(mapped_edge = rownames(tree$mapped.edge)[which(temp_map==1)])
    to_return$a <- a
    to_return$b <- b
    
    return(to_return)
  }
  
  for(i in 1:nrow(non_zero)){
    state1 <- states[non_zero$row[i]]
    state2 <- states[non_zero$col[i]]
    res[[i]] <- aa(tree, state1, state2)
  }
  
  results <- do.call(rbind, res)
  
  return(results)
}

res <- bb(simmap_best[[1]])

head(res)

res

# do it for 1000 trees
res <- list()

pb <- progress_bar$new(total = length(simmap_best))

for(i in 1:length(simmap_best)){
  pb$tick()
  res[[i]] <- bb(simmap_best[[i]])
}

res <- bind_rows(res, .id = "index")

# save out this object to read back in
saveRDS(res, here('data/sequencing_rpoB/processed/transition_rates/time_of_transition.rds'))

