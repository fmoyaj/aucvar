test_that("multiplication works", {
  expect_equal(
    {
      set.seed(100)
      library(VGAM)
      N <- 500
      true.rho <- 0.7
      data.mat <- rbinorm(N, cov12 =  true.rho)  # Bivariate normal
      x <- data.mat[, 1]
      y <- data.mat[, 2]
      u_stat_value <- kendall.tau(x,y,exact=TRUE) # Exact value of Kendall Tau statistic
      # Phi function for the Kendall Tau statistic
      kernel <- function (s1, s2){
      result <- 0
      ind <- 0
      if (s1[1] < s2[1] && s1[2] < s2[2]){
      ind <- 1
      } else if (s1[1] > s2[1] && s1[2] > s2[2]) {
      ind <- 1
      }
      return (2*(ind)-1)
      }
      varUn(x, y, 2, 2, kernel, 10^3, u_stat_value)
    }, 0.004)
})
