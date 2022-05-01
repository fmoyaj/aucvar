#' Asymptotic confidence interval of AUC based on an unbiased variance estimator
#' and the asymptotic normality
#'
#' @param p_pred A vector of the predicted probabilities for the observations
#' in the data set
#' @param label_true A vector of the true labels in the dataset, coded as 1
#' (positive) and 0 (negative)
#' @param conf_level The confidence level required.The default is 0.95.
#' @param B An integer indicating the desired number of bootstrap samples to
#' calculate the variance of the AUC. If B is set to Inf or if it is ommited,
#' the exact number of possible partitions will be calculated.
#'
#' @return A matrix (or vector) with columns giving lower and upper confidence limits
#' for each parameter. These will be labeled as (1-level)/2 and 1 - (1-level)/2 in %
#' (2.5% and 97.5% by default).

#' @export
#'
#' @examples
#' library(aucvar)
#' mydata <- na.omit(breastcancer) # Omit NA values
#' optimal_model <- glm(Class~`Clump Thickness`+`Uniformity of Cell Shape`+
#' `Bare Nuclei` + `Bland Chromatin`, family=binomial(link="logit"), data=mydata)
#' prob <- predict(optimal_model, type="response")
#' labels <- mydata$Class
#' CI_AUC(prob, labels)
CI_AUC <- function(p_pred, label_true, conf_level = 0.95, B = Inf)
{
  # Checking arguments
  # conf_level must be between 0 and 1
  if (conf_level < 0 || conf_level > 1)
  {
    base::stop("conf_level must be between 0 and 1")
  }

  # B must be a number or Inf
  if (B != Inf & B%%1 != 0)
  {
    base::stop("B must be equal to an integer number, Inf or not specified.")
  }

  auc.value <- aucvar::auc(p_pred, label_true)
  se <- base::sqrt(aucvar::varAUC(p_pred, label_true, B))
  margin.error <- stats::qnorm(conf_level + (1-(conf_level))/2) * se
  result <- base::matrix(c(auc.value - margin.error, auc.value + margin.error),
                         nrow = 1, dimnames = list(c(NULL),
                                                   c(paste(toString((1-conf_level)*100/2), '%'),
                                                     paste(toString(100-(1-conf_level)*100/2), '%'))))
  return (result)
}
