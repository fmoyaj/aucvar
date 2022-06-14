# Helper function
# Convert a vector of labels into positive (1) and negative classes (0)
makeBinaryLabels <- function(label_data)
{
  labels <- base::as.factor(label_data)
  if (base::all(levels(labels) != c(0,1)))
  { # Rename levels of factor if necessary
    base::levels(labels) <- c(0,1)
  }
  return(labels)
}


exactVarAUC <- function(label_true, p_pred, n1, n2)
{
  pos.label <-  utils::combn(n1,2)
  neg.label <-  utils::combn(n2,2)

  phat.pos <- p_pred[label_true == 1]
  phat.neg <- p_pred[label_true == 0]

  Q0 <- 0
  Q2 <- aucvar::auc(p_pred, label_true)^2 # Possible source of errors.Check whether auc works with labels

  N0 <- 2*base::dim(pos.label)[1] * base::dim(neg.label)[1]

  for(i in 1:base::dim(pos.label)[1])
  {
    for(j in 1:base::dim(neg.label)[1])
    {

      Q0 <- Q0 + (1/N0) * base::as.numeric(phat.pos[pos.label[i,1]] >
                                             phat.neg[neg.label[j,1]]) *
        base::as.numeric(phat.pos[pos.label[i,2]] > phat.neg[neg.label[j,2]]) +
        (1/N0) * base::as.numeric(phat.pos[pos.label[i,1]] > phat.neg[neg.label[j,2]]) *
        base::as.numeric(phat.pos[pos.label[i,2]] > phat.neg[neg.label[j,1]])
    }
  }

  return (Q2 - Q0)
}


#' Calculate an unbiased variance estimator of AUC
#'
#' @param p_pred A vector of the predicted probabilities for the observations in the data set
#' @param label_true A vector of the true labels in the data set, coded as 1 (positive) and 0 (negative)
#' @param B The number of random partitions to use in the partition-resampling scheme. If B is
#' set to Inf or it is ommited, the exact unbiased variance estimation formula without the partition
#' resampling scheme is realized.
#'
#' @return The value of the AUC variance estimator
#' @export
#'
#' @examples
#' library(aucvar)
#' my_data <- na.omit(breastcancer) # Omit NA values
#' my_model <- glm(Class~`Clump Thickness`+ `Uniformity of Cell Size` +
#' `Uniformity of Cell Shape`+ `Marginal Adhesion` + `Single Epithelial Cell Size` +
#' `Bare Nuclei` + `Bland Chromatin` + `Normal Nucleoli` + `Mitoses`,
#' family=binomial(link="logit"), data=my_data)
#' predictions <- predict(my_model, type="response")
#' varAUC(predictions, my_data$Class, 10^3)
#'
#' @references
#' \cite{Q. Wang and A. Guo (2020). An efficient variance estimator of AUC with applications to
#' binary classification. Statistics in Medicine 39 (28): 4281-4300. DOI: 10.1002/sim.8725.}
#'
varAUC <- function(p_pred,
                   label_true,
                   B = Inf){

  # Checking arguments
  # label_true must be a vector or a factor
  if (!base::is.vector(label_true) && !base::is.factor(label_true))
  {
    base::stop("label_true must be a vector or a factor")
  }

  # label_true must have two levels
  if (!base::length(base::levels(base::as.factor(label_true)))){
    base::stop("label_true must have two levels")
  }

  # B must be an integer
  if (B%%1 != 0 & B != Inf)
  {
    base::stop("B must be an integer number or Inf")
  }

  # Convert label_true vector into a factor
  labels <- makeBinaryLabels(label_true)

  # Identify the minimum number of observations out of positive and negative classes
  posLabelLen <- base::length(base::subset(labels, labels==1))
  negLabelLen <- base::length(base::subset(labels, labels==0))
  m <- base::min(posLabelLen, negLabelLen)


  # If B is Inf or not present, calculate the exact number of partitions
  # Go with the original definition of the formula Q(2) - Q(0)
  if (B == Inf)
  {
    var <- exactVarAUC(labels, p_pred, posLabelLen, negLabelLen)
    return (var)
  }else
  {
    phi_bar_b <-base::array(NA,B)
    WPSS <-base::rep(NA,B)

    for (b in 1:B)
    {
      # For each sampled partition
      kernel <-base::array(NA,m)

      pos_shuf <- base::sample(base::which(labels == 1), m, replace = F)
      neg_shuf <- base::sample(base::which(labels == 0), m, replace = F)


      # Calculate the kernel function for each block of paired data
      for (pair in 1:m)
      {
        kernel[pair] <- base::as.numeric(p_pred[pos_shuf[pair]] > p_pred[neg_shuf[pair]])
      }

      phi_bar_b[b] <- base::mean(kernel)

      # Compute the within-partition sum of squares
      WPSS[b] <- (1/(m-1)) * base::sum((kernel-phi_bar_b[b])^2)

    }

    # Compute the between-partition sum of squares
    BPSS <- base::mean((phi_bar_b - mean(phi_bar_b))^2)

    # Return the partition-resampling variation of the AUC variance
    return ((1/m)*base::mean(WPSS) - BPSS)
  }
}
