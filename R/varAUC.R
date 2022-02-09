# Possible errors:
# label_true is not valid (does not have two level, has NAs)



varAUC <- function(p_pred, label_true, B){
  # Convert label_true vector into a factor
  labels <- as.factor(label_true)
  if (all(levels(labels) != c(0,1))){ # Rename levels of factor if necessary
    levels(labels) <- c(0,1)
  }

  # Minimum number of observations out of positive and negative classes
  m <- min(length(subset(labels, labels==1)), length(subset(labels, labels==0)))

  phi_bar_b <-array(NA,B)
  WPSS <-rep(NA,B)

  for (b in 1:B){
    # For each sampled partition
    kernel <-array(NA,m)

    pos_shuf <-sample(which(labels == 1), m, replace = F)
    neg_shuf <-sample(which(labels == 0), m, replace = F)


    # Calculate the kernel function for each block of paired data
    for (pair in 1:m){
      kernel[pair] <- as.numeric(p_pred[pos_shuf[pair]] > p_pred[neg_shuf[pair]])
    }

    phi_bar_b[b] <- mean(kernel)

    # Compute the within-partition sum of squares
    WPSS[b] <- (1/(m-1)) * sum(kernel-phi_bar_b[b])^2

  }

  # Compute the between-partition sum of squares
  BPSS <- (1/B) * mean((phi_bar_b - mean(phi_bar_b))^2)

  # Return the partition-resampling variation of the AUC variance
  return (1/m)*mean(WPSS) - BPSS
}
