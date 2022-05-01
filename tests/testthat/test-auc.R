test_that("AUC calculation with optimal model for breast cancer data", {
  expect_equal({
    library(aucvar)
    mydata <- na.omit(breastcancer)
    optimal_model <- glm(Class~`Clump Thickness`+`Uniformity of Cell Shape`+
      `Bare Nuclei` + `Bland Chromatin`, family=binomial(link="logit"), data=mydata)
    prob <- predict(optimal_model, type="response")
    labels <- mydata$Class
    auc(prob, labels)
  }, 0.99520807)
})

test_that("AUC calculation with model with uniformity of cell shape", {
  expect_equal({
    library(aucvar)
    mydata <- na.omit(breastcancer)
    optimal_model <- glm(Class~`Uniformity of Cell Shape`
                         , family=binomial(link="logit"), data=mydata)
    prob <- predict(optimal_model, type="response")
    labels <- mydata$Class
    auc(prob, labels)
  }, 0.97542783)
})
