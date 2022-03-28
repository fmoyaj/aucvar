# Possible errors:
# label_true is not valid (does not have two level, has NAs)
# data not in adequate format
# Catch errors from other functions e.g. auc

# Private function
# Convert a vector of labels into positive (1) and negative classes (0)
makeBinaryLabels <- function(label_data){
  labels <- base::as.factor(label_data)
  if (base::all(levels(labels) != c(0,1))){ # Rename levels of factor if necessary
    base::levels(labels) <- c(0,1)
  }
  return(labels)
}


exactVarAUC <- function(labels, p_pred, m, n1, n2){
  pos.label <-  utils::combn(n1,2)
  neg.label <-  utils::combn(n2,2)

  phat.pos <- p_pred[labels == 1]
  phat.neg <- p_pred[labels == 0]

  Q0 <- 0
  Q2 <- aucvar::auc(p_pred, labels)^2 # Possible source of errors.Check whether auc works with labels

  N0 <- 2*base::dim(pos.label)[1] * base::dim(neg.label)[1]

  for(I in 1:base::dim(pos.label)[1]){
    for(j in 1:base::dim(neg.label)[1]){

        Q0 <- Q0 + (1/N0) * base::as.numeric(phat.pos[pos.label[I,1]] >
                                               phat.neg[neg.label[j,1]]) *
            base::as.numeric(phat.pos[pos.label[I,2]] > phat.neg[neg.label[j,2]]) +
            (1/N0) * base::as.numeric(phat.pos[pos.label[I,1]] > phat.neg[neg.label[j,2]]) *
            base::as.numeric(phat.pos[pos.label[I,2]] > phat.neg[neg.label[j,1]])
    }
  }

  return (Q2 - Q0)
}


#' Calculate an unbiased variance estimator of AUC
#'
#' @param p_pred A vector of the predicted probabilities for the observations in the data set
#' @param label_true A vector of the true labels in the dataset, coded as 1 (positive) and 0 (negative)
#' @param B The number of random partitions to use in the partition-resampling scheme. If B is
#' set to Inf or it is ommited, the exact number of possible partitions will be calculated.
#'
#' @return The value of the AUC variance estimator
#' @export
#'
#' @examples
#' library(aucvar)
#' data <- na.omit(breastcancer) # Omit NA values
#' optimal_model <- glm(Class~`Clump Thickness`+`Uniformity of Cell Shape`+
#' `Bare Nuclei` + `Bland Chromatin`, family=binomial(link="logit"), data=data)
#' predictions <- predict(optimal_model, type="response")
#' varAUC(predictions, data$Class, 10)
varAUC <- function(p_pred,
                   label_true,
                   B = Inf){

  # Convert label_true vector into a factor
  labels <- makeBinaryLabels(label_true)

  # Identify the minimum number of observations out of positive and negative classes
  posLabelLen <- base::length(base::subset(labels, labels==1))
  negLabelLen <- base::length(base::subset(labels, labels==0))
  m <- base::min(posLabelLen, negLabelLen)

  # If B is Inf or not present, calculate the exact number of partitions
  # Go with the original definition of the formula Q(2) - Q(0)
  if (B == Inf){
      var <- exactVarAUC(labels, p_pred, m, posLabelLen, negLabelLen)
      return (var)
  }

  phi_bar_b <-base::array(NA,B)
  WPSS <-base::rep(NA,B)

  for (b in 1:B){
    # For each sampled partition
    kernel <-base::array(NA,m)

    pos_shuf <- base::sample(base::which(labels == 1), m, replace = F)
    neg_shuf <- base::sample(base::which(labels == 0), m, replace = F)


    # Calculate the kernel function for each block of paired data
    for (pair in 1:m){
      kernel[pair] <- base::as.numeric(p_pred[pos_shuf[pair]] > p_pred[neg_shuf[pair]])
    }

    phi_bar_b[b] <- base::mean(kernel)

    # Compute the within-partition sum of squares
    WPSS[b] <- (1/(m-1)) * base::sum(kernel-phi_bar_b[b])^2

  }

  # Compute the between-partition sum of squares
  BPSS <- (1/B) * base::mean((phi_bar_b - mean(phi_bar_b))^2)

  # Return the partition-resampling variation of the AUC variance
  return (1/m)*base::mean(WPSS) - BPSS
}
