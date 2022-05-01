
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aucvar

<!-- badges: start -->
<!-- badges: end -->

The goal of aucvar is to calculate the variance of AUC through different
estimators, including an unbiased variance estimator, bootstrap variance
estimator, and jackknife variance estimator.

## Installation

You can install the development version of aucvar from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fmoyaj/aucvar")
```

## Example

This is a basic example which shows you how to calculate the variance of
the AUC using an unbiased estimator:

``` r
library(aucvar)
mydata <- na.omit(breastcancer) # Omit NA values
optimal_model <- glm(Class~`Clump Thickness`+`Uniformity of Cell Shape`+
      `Bare Nuclei` + `Bland Chromatin`, family=binomial(link="logit"), data=mydata)
predictions <- predict(optimal_model, type="response")
varAUC(predictions, mydata$Class, 10^3)
#> [1] 2.875184e-06
```

You can also calculate the AUC for a logistic regression model

``` r
library(aucvar)
mydata <- na.omit(breastcancer) # Omit NA values
optimal_model <- glm(Class~`Uniformity of Cell Shape`
                         , family=binomial(link="logit"), data=mydata)
prob <- predict(optimal_model, type="response")
labels <- mydata$Class
auc(prob, labels)
#> [1] 0.9754278
```
