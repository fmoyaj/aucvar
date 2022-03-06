# fomula look at ?formula from stats package
# allow to change link? - YES. TODO: Add error
# data needs to be dataset without the response
# formula is not in the right format
# takes too much time to exectute
# how many decimal points for SD ?
# Figure out deal with helper functions to not have them in two separate files
# TODO:example


# Private function
# Convert a vector of labels into positive (1) and negative classes (0)
makeBinaryLabels <- function(label_data){
  labels <- base::as.factor(label_data)
  if (base::all(levels(labels) != c(0,1))){ # Rename levels of factor if necessary
    base::levels(labels) <- c(0,1)
  }
  return(labels)
}




#' Calculate the bootstrap variance estimator for AUC
#'
#' @param formula A string with an expression of the form `y ~ model` that represents
#' the binary classification model. It may include operators as +, ^ and :
#' @param label_true A vector of the true labels in the data set, coded as 1 (positive) and 0 (negative)
#' @param data A data frame, list or environment containing the variables in the model
#' except for the response variable. It can also be an object coercible by as.data.frame
#' to a data frame.
#' @param B An integer indicating the desired number of bootstrap samples
#' @param link A string specifying the model link function for glm function
#' used to fit the binomial model. Possible links include `logit`., `probit`, `cauchit`,
#' (corresponding to logistic, normal and Cauchy CDFs respectively) `log` and
#' `cloglog` (complementary log-log). The default is `logit`.
#' @return The value of the bootstrap variance estimator for the AUC.
#' @export
#'
#' @examples
#' 'TODO'
var_boot <- function(formula, label_true, data, B, link = "logit"){

  # Convert label_true vector into a factor
  labels <- makeBinaryLabels(label_true)

  # Convert string formula into a formula object from the stats package
  formula.obj <- stats::formula(formula)

  posLabelLen <- base::length(base::subset(labels, labels==1))
  negLabelLen <- base::length(base::subset(labels, labels==0))

  boot.AUC <-base::array(NA, B)

  N <- base::nrow(data) # Size of data, which is also the size of the sample

  # For each sampled partition
  for (b in 1:B){
    boot.sample <-data[base::sample(N, N, replace = T),] # Draw sample of size N with replacement

    # If the sample is all O's or all 1's, assume that sample is improbable
    if (base::unique(boot.sample) == 1){
      next
    }

    # Update entries in boot.AUC array with calculated AUC values
    # Choose the link, add options: ?????
    fit <-stats::glm(formula.obj, family = stats::binomial(link = link), data = boot.sample)
    boot.AUC[b] <-aucvar::auc(stats::predict.glm(fit, newdata = boot.sample, type = "response"))
  }

  return(stats::var(boot.AUC))
}
