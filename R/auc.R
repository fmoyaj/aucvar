# Add error if function stops/returns null
# If there are more than two distinct label symbols, execution stops with an error message. -> propagate exception

#' Calculate the area under the ROC curve
#'
#' @param p_pred A vector of the predicted probabilities for the observations
#' in the data set.
#' @param label_true A vector containing the true class labels in the dataset.
#'
#' @return Area under the ROC curve, which is equal to the Wilcoxon-Mann-Whitney
#'test statistic and also the probability that the classifier will score are
#'randomly drawn positive sample higher than a randomly drawn negative sample.
#'
#' @export
#'
#' @examples
#' library(aucvar)
#' data <- na.omit(breastcancer) # Omit NA values
#' optimal_model <- glm(Class~`Clump Thickness`+`Uniformity of Cell Shape`+
#' `Bare Nuclei` + `Bland Chromatin`, family=binomial(link="logit"), data=data)
#' prob <- predict(optimal_model, type="response")
#' labels <- data$Class
#' auc(prob, labels)
auc <- function(p_pred, label_true){
  pred <-ROCR::prediction(p_pred, label_true)
  perf <- ROCR::performance(pred, measure = "auc", fpr.stop=1) # Stop if FPR is 1
  return(base::as.numeric(perf@y.values))
}
