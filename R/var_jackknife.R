# Check whether N > d, if not, report an error
# Add recommendations for d in the recommendation
# Special case:if d = 1,that is the default. Total number of n possibilities. Add delete-1 version.
# How to modify code for delete 1 version? Delete each observation, one at a time. Fit model n times???
# Add warning for delete 1 -> function needs to be smooth enough




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


#' Calculate the jackknife variance estimator for AUC
#'
#' @param formula_string A string with an expression of the form y ~ model that represents
#' the binary classification model. It may include operators as +, ^ and :
#' @param label_true A vector of the true labels in the data set, coded as 1 (positive) and 0 (negative)
#' @param data A data frame, list or environment containing the variables in the model
#' except for the response variable. It can also be an object coercible by as.data.frame
#' to a data frame.
#' @param B An integer indicating the desired number of jackknife samples.
#' @param d Number of data entries to remove to generate the jackknife samples.
#' If d is equal to 1, the delete-one version of the jackknife variance estimator
#' will be used and B does not need to be specified.
#' @param link A string specifying the model link function for glm function
#' used to fit the binomial model. Possible links include `logit`., `probit`, `cauchit`,
#' (corresponding to logistic, normal and Cauchy CDFs respectively) `log` and
#' `cloglog` (complementary log-log). The default is `logit`.
#'
#' @return The value of the jackknife variance estimator for the AUC
#' @export
#'
#' @examples
#' library(aucvar)
#' data <- na.omit(breastcancer) # Omit NA values
#' model_formula <- "Class~`Clump Thickness`+`Uniformity of Cell Shape`+
#' `Bare Nuclei` + `Bland Chromatin`" # Use quotes inside double quotes since
#' # dataset variable names have spaces
#' var_jackknife(model_formula, data$Class, data, B = 10^3, d = 20)
var_jackknife <- function(formula_string, label_true, data, B = Inf, d, link = "logit")
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

  N <- base::nrow(data) # Size of data, which is also the size of the sample

  # Calculating the total number of partitions if B is Inf or not provided -
  if (B == Inf)
  {
    B <- base::choose(N, d)
  }

  jack.AUC <-base::array(NA, B)

  for (b in 1:B){
    index <- base::sample(N, N-d, replace = F)
    jack.sample <-data[index,] # Draw sample of size N-d with no replacement

    # If the sample is all O's or all 1's, assume that sample is improbable
    if (all(labels[index] == 0) || all(labels[index] == 1))
    {
      next
    }

    # Update entries in jack.AUC array with calculated AUC values
    fit <-stats::glm(formula.obj, family = stats::binomial(link = link), data = jack.sample)
    jack.AUC[b] <- aucvar::auc(stats::predict.glm(fit, newdata = jack.sample, type = "response"),
                               labels[index])


  }

  return (((N-d)/(d*B))* base::sum((jack.AUC-mean(jack.AUC))^2))


}
