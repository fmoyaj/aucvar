
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aucvar

<!-- badges: start -->
<!-- badges: end -->

AUC, area under an ROC curve, is one of the most commonly used measures
to evaluate the performance of a binary classifier across various
discrimination thresholds between 0 and 1 (Bradley 1997). The developed
R package “aucvar” realizes several variance estimation methods for AUC
based on its two-sample U-statistic expression. It complements existing
R packages that mainly focus on the calculation of AUC or visualization
of the ROC curve alone. In particular, it aims at facilitating
statistical inference of AUC in practical applications. In addition to
computing the unbiased variance estimator of AUC proposed by Wang and
Guo (2020), “aucvar” also offers tools that can be used to realize the
naive nonparametric bootstrap (Efron 1979) and the delete-d jackknife
(Efron and Stein 1981) variance estimators. Moreover, there is a
ready-to-use function to construct an asymptotic confidence interval of
AUC. As an extension to the AUC-specific application, the package also
includes a function that computes the unbiased variance estimator of a
general two-sample U-statistic.

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
my_data <- na.omit(breastcancer) # Omit NA values
full_model <- glm(Class~`Clump Thickness`+ `Uniformity of Cell Size` +
`Uniformity of Cell Shape`+ `Marginal Adhesion` + `Single Epithelial Cell Size` +
`Bare Nuclei` + `Bland Chromatin` + `Normal Nucleoli` + `Mitoses`,
family=binomial(link="logit"), data=my_data)
predictions <- predict(full_model, type="response")
varAUC(predictions, my_data$Class, 10^3)
#> [1] 1.726719e-06
```

You can also calculate the AUC for a logistic regression model

``` r
library(aucvar)
my_data <- na.omit(breastcancer) # Omit NA values
full_model <- glm(Class~`Clump Thickness`+ `Uniformity of Cell Size` +
`Uniformity of Cell Shape`+ `Marginal Adhesion` + `Single Epithelial Cell Size` +
`Bare Nuclei` + `Bland Chromatin` + `Normal Nucleoli` + `Mitoses`,
family=binomial(link="logit"), data=my_data)
prob <- predict(full_model, type="response")
labels <- my_data$Class
auc(prob, labels)
#> [1] 0.9963248
```
