test_that("Unbiased variance calculation with optimal model for breast cancer data", {
  expect_equal(
    {
      set.seed(100)
      library(aucvar)
      mydata <- na.omit(breastcancer) # Omit NA values
      optimal_model <- glm(Class~`Clump Thickness`+`Uniformity of Cell Shape`+
      `Bare Nuclei` + `Bland Chromatin`, family=binomial(link="logit"), data=mydata)
      predictions <- predict(optimal_model, type="response")
      varAUC(predictions, mydata$Class, 10^3)
    }, 1.41465743e-06)
})
