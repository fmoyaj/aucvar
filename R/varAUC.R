# Possible errors:
# label_true is not valid (does not have two level, has NAs)



#' Calculate an unbiased variance estimator of AUC
#'
#' @param p_pred A vector of the predicted probabilities for the observations in the data set
#' @param label_true A vector of the true labels in the dataset, coded as 1 (positive) and 0 (negative)
#' @param B The number of random partitions to use in the partition-resampling scheme
#'
#' @return A float indicating the value of the AUC variance estimator
#' @export
#'
#' @examples
#' varAUC(c(0.50,0.36,0.68,0.10,0.54,0.36,0.88,0.95,0.32,0.45,0.12,0.97,0.66,
#' 0.22,0.45,0.19,0.87,0.44,0.32,0.11,0.94,0.43,0.55,0.32,0.67,0.67,0.23),
#' c(1,0,1,1,0,0,1,0,1,1,1,0,1,1,0,1,0,1,0,1,1,1,1,1,1,1,0), 10)
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
