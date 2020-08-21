# Nothing here is exported

# partition ---------------------------------------------------------------

# Creates the time interval partition
partition <- function(t, int_len) {
  
  last <- ifelse(max(t) %% int_len == 0, max(t), max(t) + int_len)
  seq(from = 0, to = last, by = int_len)
  
}


# partition_location ------------------------------------------------------

# Returns the interval of the partition in which t is located
partition_location <- function(t, partition) {
  
  lt <- length(t)
  t_loc <- rep(0, times = lt)
  
  for (i in 1:lt) {
    j <- 1
    while(t[i] > partition[j + 1] & j < length(partition)) {
      j <- j + 1
    }
    t_loc[i] <- j
  }
  
  return(t_loc)
  
}


# partition_count ---------------------------------------------------------

# Count how many times are at each interval
partition_count <- function(t, partition) {
  
  t_loc <- partition_location(t, partition)
  counts <- rep(0, times = (length(partition) - 1))
  
  for (i in 1:length(counts)) {
    counts[i] <- sum(t_loc == i)
  }
  
  return(counts)
  
}


# cum_h -------------------------------------------------------------------

# Compute the cumulative hazard H(t)
cum_h <- function(t, partition, lambda) {
  
  cumhaz <- vector(mode = "double", length = length(t))
  int_len <- partition[2] - partition[1]
  index <- partition_location(t, partition)
  
  for (i in 1:length(t)) {
    loc <- index[i]
    if (loc > length(lambda)) {
      loc <- length(lambda)
    }
    cum <- 0
    if (loc > 1) {
      cum <- sum(lambda[1:(loc - 1)]) * int_len
    }
    cumhaz[i] <- cum + (t[i] - partition[loc]) * lambda[loc]
  }
  
  return(cumhaz)
  
}
