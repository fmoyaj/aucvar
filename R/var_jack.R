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

# Delete-one version of jackknife
deleteOne <- function(formula.obj, labels, data, link) {

  N <- base::nrow(data) # Size of data, which is also the size of the sample
  jack.AUC <-base::array(NA, N)

  for (i in 1:N)
    {

    jack.sample <-data[-i,] # Take out ith observation

    # Update entries in jack.AUC array with calculated AUC values
    fit <-stats::glm(formula.obj, family = stats::binomial(link = link), data = jack.sample)
    jack.AUC[i] <- aucvar::auc(stats::predict.glm(fit, newdata = jack.sample, type = "response"),
                               labels[-i])


  }

  return (((N-1)/(N))* base::sum((jack.AUC-mean(jack.AUC))^2))
}


#' Calculate the delete-d jackknife variance estimator for AUC
#'
#' @param formula_string A string with an expression of the form y ~ model that represents
#' the binary classification model. It may include operators as +, ^ and :
#' @param label_true A vector of the true labels in the data set, coded as 1 (positive) and 0 (negative)
#' @param data A data frame, list or environment containing the variables in the model.
#' It can also be an object coercible by as.data.frame
#' to a data frame.
#' @param B An integer indicating the desired number of jackknife samples.
#' @param d Number of data entries to remove to generate the jackknife samples.
#' If d is equal to 1, the delete-one version of the jackknife variance estimator
#' will be used in which case B does not need to be specified.
#' @param link A string specifying the model link function for glm function
#' used to fit the binomial model. Possible links include `logit`., `probit`, `cauchit`,
#' (corresponding to logistic, normal and Cauchy CDFs respectively) `log` and
#' `cloglog` (complementary log-log). The default is `logit`.
#'
#' @return The value of the delete-d jackknife variance estimator for the AUC
#' @export
#'
#' @examples
#' library(aucvar)
#' my_data <- na.omit(breastcancer) # Omit NA values
#' model_formula <- "Class~`Clump Thickness`+ `Uniformity of Cell Size`+`Uniformity of Cell Shape`+
#' `Marginal Adhesion` + `Single Epithelial Cell Size` + `Bare Nuclei` +
#' `Bland Chromatin` + `Normal Nucleoli` + `Mitoses`"
#' # Use quotes inside double quotes since data set variable names have spaces
#' var_jack(model_formula, my_data$Class, my_data, B = 10^3, d = 20)
#'
#' @references
#' \cite{B. Efron (1979). Bootstrap methods: another look at the jackknife The Annals of
#' Statistics 7: 1-26.}
#'
var_jack <- function(formula_string, label_true, data, B = Inf, d, link = "logit")
{
  # Check arguments
  # formula_string must be a string
  if (base::typeof(formula_string) != "character")
  {
    base::stop("formula_string must be a string")
  }

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
  if (B%%1 != 0)
  {
    base::stop("B must be an integer number")
  }

  # Links must be in the list of accepted links for the binomial family
  if (!(link %in% base::c("logit", "probit", "cauchit", "log", "cloglog")))
  {
    base::stop("Link provided is not valid. Link must be logit, probit,
         cauchit, log, cloglog")
  }

  N <- base::nrow(data) # Size of data, which is also the size of the sample

  # N must be greater than d
  if (N < d)
  {
    base::stop("d must be smaller than the length of the data set provided")
  }

  # Convert label_true vector into a factor
  labels <- makeBinaryLabels(label_true)

  # Convert string formula into a formula object from the stats package
  formula.obj <- base::tryCatch({stats::formula(formula_string)}, error = function(){
    base::stop("Could not convert given formula_string into a formula object.")
  })


  # If d is equal to 1, do delete-one version of jackknife calculation and return
  if (d == 1)
  {
    return (deleteOne(formula.obj, labels, data, link))
  }

  posLabelLen <- base::length(base::subset(labels, labels==1))
  negLabelLen <- base::length(base::subset(labels, labels==0))

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
