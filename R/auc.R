#' Calculate the area under the ROC curve
#'
#' @param p_pred A vector of the predicted probabilities for the observations
#' in the data set.
#' @param label_true A vector containing the true class labels in the dataset.
#'
#' @return Area under the ROC curve, which is equal to the Mann-Whitney
#'two-sample U-statistic. It is also the estimated probability that the binary
#'classifier will score a randomly drawn positive sample higher than a
#'randomly drawn negative sample.
#'
#' @export
#'
#' @examples
#' library(aucvar)
#' mydata <- na.omit(breastcancer) # Omit NA values
#' optimal_model <- glm(Class~`Clump Thickness`+`Uniformity of Cell Shape`+
#' `Bare Nuclei` + `Bland Chromatin`, family=binomial(link="logit"), data=mydata)
#' prob <- predict(optimal_model, type="response")
#' labels <- mydata$Class
#' auc(prob, labels)
auc <- function(p_pred, label_true){

  # Check arguments
  if((!is.vector(p_pred) || !is.vector(label_true)) && !is.factor(label_true))
  {
    base::stop("p_pred and label_true must be vectors or a factor")
  }

  # label_true must have the same length as p_pred
  pred_length <- base::length(p_pred)
  label_length <- base::length(label_true)
  if(pred_length != label_length){
    base::stop("p_pred and label_true must have the same length")
  }

  # label_true must have two levels
  if (!base::length(base::levels(base::as.factor(label_true)))){
    base::stop("label_true must have two levels")
  }

  pred <-ROCR::prediction(p_pred, label_true)
  perf <- ROCR::performance(pred, measure = "auc", fpr.stop=1) # Stop if FPR is 1
  return(base::as.numeric(perf@y.values))
}
