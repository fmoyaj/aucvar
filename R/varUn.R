# Add exact definition when B = Inf

varUn <- function(sample1, sample2, k1, k2, phi ,B = Inf){

# WARNING: k1, k2 cannot be greater than half of the corresponding sample
  m1 <- base::nrow(sample1)
  m2 <- base::nrow(sample2)

  # Maximum number of partitions we can get
  m <- base::min(base::floor(m1/k1), base::floor(m2/k2))

  for (i in 1:B){
    # For each sampled partition
    kernel <-base::array(NA,m)

    id1.sub <- base::array(NA,m)
    id2.sub <- base::array(NA,m)

    # Shuffle k indexes for m blocks
    for (block in 1:m){
      id1.sub[block] <- base::sample(m1, k1, replace = F)
      id2.sub[block] <- base::sample(m2, k2, replace = F)
    }

    # Calculate the kernel function for each block of paired data
    for (block in 1:m){
      kernel[block] <- phi(sample1[id1.sub,],sample2[id2.sub,])
    }

    phi_bar_b[b] <- base::mean(kernel)

    # Compute the within-partition sum of squares
    WPSS[b] <- (1/(m-1)) * base::sum(kernel-phi_bar_b[b])^2
  }

  # Return the partition-resampling variation of the AUC variance
  return (1/m)*base::mean(WPSS) - BPSS
}
