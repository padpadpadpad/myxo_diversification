# lets see if we can play around with the stochastic character mapping

library(phytools)

data(anoletree)
anoletree
tree <- anoletree

# try and write an example to grab the mapped edges of all transitions
# grab out code from countSimmap()

# calculate total number of transitions
n <- sum(sapply(tree$maps, length)) - nrow(tree$edge)

# set states
if (is.null(states)) 
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

results
