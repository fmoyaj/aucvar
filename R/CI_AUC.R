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


#' Asymptotic confidence interval of AUC based on a variance estimator
#' and the asymptotic normality
#'
#' @param formula_string A string with an expression of the form `y ~ model` that represents
#' the binary classification model. It may include operators as +, ^ and :
#' @param link A string specifying the model link function for glm function
#' used to fit the binomial model. Possible links are `logit`., `probit`, `cauchit`,
#' (corresponding to logistic, normal and Cauchy CDFs respectively) `log` and
#' `cloglog` (complementary log-log). The default is `logit`.
#' @param data A data frame, list or environment containing the variables in the model
#' except for the response variable. It can also be an object coercible by as.data.frame
#' to a data frame.
#' @param label_true A vector of the true labels in the data set, coded as 1
#' (positive) and 0 (negative)
#' @param conf_level The confidence level required.The default is 0.95.
#' @param method The method to use to compute the variance estimator. The possible
#' methods are "unbiased" for the unbiased variance estimator of the AUC devised by
#' Wang and Guo, "jackknife" for the jackknife variance estimator, and "bootstrap"
#' for the bootstrap variance estimator of the AUC.
#' @param B if the method chosen is "unbiased", B is an integer indicating the
#' desired number of bootstrap samples to calculate the variance of the AUC.
#' If the method chosen is "unbiased" or "jackknife" and B is set to Inf or if it is omitted,
#' the exact number of possible partitions will be calculated.
#' @param d Number of data entries to remove to generate the jackknife samples for the
#' jackknife variance estimator. This argument is only required when the method chosen
#' is "jackknife". If d is equal to 1, the delete-one version of the jackknife
#' variance estimator will be used in which case B does not need to be specified.
#'
#' @return A matrix (or vector) with columns giving lower and upper confidence limits
#' for each confidence level. These will be labeled as (1-level)/2 and 1 - (1-level)/2 in
#' percentage. (2.5% and 97.5% by default).
#'
#' @export
#'
#' @examples
#' library(aucvar)
#' my_data <- na.omit(breastcancer) # Omit NA values
#' model_formula <- "Class~`Clump Thickness`+ `Uniformity of Cell Size`+`Uniformity of Cell Shape`+
#' `Marginal Adhesion` + `Single Epithelial Cell Size` + `Bare Nuclei` +
#' `Bland Chromatin` + `Normal Nucleoli` + `Mitoses`"
#' # Use quotes inside double quotes since data set variable names have spaces
#' labels <- my_data$Class
#' CI_AUC(model_formula, "logit", my_data, labels, 0.95, "unbiased", B = 10^3)
CI_AUC <- function(formula_string, link = "logit", data, label_true, conf_level = 0.95,
                   method = "unbiased",  B = Inf, d = Inf)
{
  # Checking arguments
  # conf_level must be between 0 and 1
  if (conf_level < 0 || conf_level > 1)
  {
    base::stop("conf_level must be between 0 and 1")
  }

  # method is valid
  if(!method %in% base::c("unbiased", "bootstrap", "jackknife"))
  {
    base::stop('method must be "unbiased", "bootstrap", "jackknife"')
  }

  # Links must be in the list of accepted links for the binomial family
  if (!(link %in% base::c("logit", "probit", "cauchit", "log", "cloglog")))
  {
    base::stop("Link provided is not valid. Link must be logit, probit,
         cauchit, log, cloglog")
  }

  # B can only be Inf if method is unbiased
  if(B == Inf & method %in% c("bootstrap"))
  {
    base::stop('B cannot be Inf if method is "bootstrap"')
  }

  # d must be specified if method is jackknife
  if(d == Inf & method %in% c("jackknife"))
  {
    base::stop('d must have a value of 1 or higher if method is jackknife')
  }

  # B must be a number or Inf
  if (B != Inf & B%%1 != 0)
  {
    base::stop("B must be equal to an integer number, Inf or not specified.")
  }

  # Create formula based on string provided
  formula.obj <- base::tryCatch({stats::formula(formula_string)}, error = function(){
    base::stop("Could not convert given formula_string into a formula object.")
  })


  # Convert labels
  labels <- makeBinaryLabels(label_true)
  fit <-stats::glm(formula.obj, family = stats::binomial(link = link), data = data)
  p_pred <- stats::predict.glm(fit, newdata = data, type = "response")

  # Modifying se for different methods
  if (method == "unbiased"){
    se <- base::sqrt(aucvar::varAUC(p_pred, label_true, B)) # method is "unbiased"
  }

  else {
    if (method == "bootstrap"){
      se <- base::sqrt(aucvar::var_boot(formula_string, label_true, data, B, link))
    }
    else {
      se <- base::sqrt(aucvar::var_jack(formula_string, label_true, data, B, d))
    }

  }

  auc.value <- aucvar::auc(p_pred, label_true)
  margin.error <- stats::qnorm(conf_level + (1-(conf_level))/2) * se
  result <- base::matrix(c(auc.value - margin.error, auc.value + margin.error),
                         nrow = 1, dimnames = list(c(NULL),
                                                   c(paste(toString((1-conf_level)*100/2), '%'),
                                                     paste(toString(100-(1-conf_level)*100/2), '%'))))
  return (result)
}
