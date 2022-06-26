test_that("Data set has 10 columns", {
  expect_equal(ncol(breastcancer), 10)
})

test_that("Data set has 699 rows", {
  expect_equal(nrow(breastcancer), 699)
})
