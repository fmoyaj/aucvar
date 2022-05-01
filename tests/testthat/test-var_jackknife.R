test_that("Jackknife variance calculation with optimal model for breast cancer data", {
  expect_equal(
    {
      set.seed(100)
      library(aucvar)
      mydata <- na.omit(breastcancer) # Omit NA values
      model_formula <- "Class~`Clump Thickness`+`Uniformity of Cell Shape`+
    `Bare Nuclei` + `Bland Chromatin`"
      var_jackknife(model_formula, mydata$Class, mydata, B = 10^3, d = 20)
    }, 2.58345346e-06)
})
