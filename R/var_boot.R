# Convert a vector of labels into a factor with positive (1) and negative classes (0)
makeBinaryLabels <- function(label_data)
{
  labels <- base::as.factor(label_data)

  if (base::all(levels(labels) != c(0,1)))
  { # Rename levels of factor if necessary
    base::levels(labels) <- c(0,1)
  }
  return(labels)
}




#' Calculate the bootstrap variance estimator for AUC
#'
#' @param formula_string A string with an expression of the form `y ~ model` that represents
#' the binary classification model. It may include operators as +, ^ and :
#' @param label_true A vector of the true labels in the data set, coded as 1 (positive) and 0 (negative)
#' @param data A data frame, list or environment containing the variables in the model
#' except for the response variable. It can also be an object coercible by as.data.frame
#' to a data frame.
#' @param B An integer indicating the desired number of bootstrap samples
#' @param link A string specifying the model link function for glm function
#' used to fit the binomial model. Possible links are `logit`., `probit`, `cauchit`,
#' (corresponding to logistic, normal and Cauchy CDFs respectively) `log` and
#' `cloglog` (complementary log-log). The default is `logit`.
#' @return The value of the bootstrap variance estimator for the AUC.
#' @export
#'
#' @examples
#' library(aucvar)
#' data <- na.omit(breastcancer) # Omit NA values
#' model_formula <- "Class~`Clump Thickness`+`Uniformity of Cell Shape`+
#' `Bare Nuclei` + `Bland Chromatin`" # Use quotes inside double quotes since
#' # dataset variable names have spaces
#' var_boot(model_formula, data$Class, data, B = 10^3)
var_boot <- function(formula_string, label_true, data, B, link = "logit")
{

  # Check arguments
  # formula_string must be a string
  if (base::typeof(formula_string) != "character")
  {
    stop("formula_string must be a string")
  }

  # label_true must be a vector or a factor
  if (!base::is.vector(label_true) && !base::is.factor(label_true))
  {
    stop("label_true must be a vector or a factor")
  }

  # label_true must have two levels
  if (!base::length(base::levels(base::as.factor(label_true)))){
    stop("label_true must have two levels")
  }

  # B must be an integer
  if (B%%1 != 0)
  {
    stop("B must be an integer number")
  }

  # Links must be in the list of accepted links for the binomial family
  if (!(link %in% c("logit", "probit", "cauchit", "log", "cloglog")))
  {
    stop("Link provided is not valid. Link must be logit, probit,
         cauchit, log, cloglog")
  }

  # Convert label_true vector into a factor
  labels <- makeBinaryLabels(label_true)

  # Convert string formula into a formula object from the stats package
  formula.obj <- stats::formula(formula_string)

  posLabelLen <- base::length(base::subset(labels, labels==1))
  negLabelLen <- base::length(base::subset(labels, labels==0))

  boot.AUC <-base::array(NA, B)

  N <- base::nrow(data) # Size of data, which is also the size of the sample

  # For each sampled partition
  for (b in 1:B)
  {
    index <- base::sample(N, N, replace = T)
    boot.sample <-data[index,] # Draw sample of size N with replacement

    # If the sample is all O's or all 1's, assume that sample is improbable
    if (all(labels[index] == 0) || all(labels[index] == 1))
    {
      next
    }

    # Update entries in boot.AUC array with calculated AUC values
    fit <-stats::glm(formula.obj, family = stats::binomial(link = link), data = boot.sample)
    boot.AUC[b] <-aucvar::auc(stats::predict.glm(fit, newdata = boot.sample, type = "response"),
                              labels[index])
  }

  return(stats::var(boot.AUC))
}
