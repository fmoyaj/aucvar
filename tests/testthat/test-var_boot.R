test_that("Bootstrap variance calculation with optimal model for breast cancer data", {
  expect_equal(
  {
    set.seed(100)
    library(aucvar)
    mydata <- na.omit(breastcancer) # Omit NA values
    model_formula <- "Class~`Clump Thickness`+`Uniformity of Cell Shape`+
    `Bare Nuclei` + `Bland Chromatin`"
    var_boot(model_formula, mydata$Class, mydata, B = 10^3)
    }, 2.68836708e-06)
})
