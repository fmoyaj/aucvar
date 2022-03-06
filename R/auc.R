# Add error if function stops/returns null
# If there are more than two distinct label symbols, execution stops with an error message. -> propagate exception

#' Calculate the area under the ROC curve
#'
#' @param p_pred A vector of the predicted probabilities for the observations
#' in the data set.
#' @param label_true A vector of the true labels in the dataset, coded as 1
#' (positive) and 0 (negative).
#'
#' @return Area under the ROC curve, which is equal to the Wilcoxon-Mann-Whitney
#'test statistic and also the probability that the classifier will score are
#'randomly drawn positive sample higher than a randomly drawn negative sample.
#'
#' @export
#'
#' @examples
#' prediction_vector <- c(0.50,0.36,0.68,0.10,0.54,0.36,0.88,0.95,0.32,0.45,0.12,
#' 0.97,0.66, 0.22,0.45,0.19,0.87,0.44,0.32,0.11,0.94,0.43,0.55,0.32,0.67,0.67,0.23)
#' label_vector <- c(1,0,1,1,0,0,1,0,1,1,1,0,1,1,0,1,0,1,0,1,1,1,1,1,1,1,0)
#' auc(prediction_vector, label_vector)
auc <- function(p_pred, label_true){
  pred <-ROCR::prediction(p_pred, label_true)
  perf <- ROCR::performance(pred, measure = "auc", fpr.stop=1) # Stop if FPR is 1
  return(base::as.numeric(perf@y.values))
}
