library(testthat)

test_that("run_fastlss returns matrix with correct dims", {
  skip_if_not_installed("hrfals")
  
  # simple matrices for testing
  n <- 100  # timepoints
  q <- 3    # nuisance regressors
  V <- 4    # voxels  
  T <- 5    # trials
  
  A <- matrix(rnorm(n * q), n, q)  # design matrix
  Y <- matrix(rnorm(n * V), n, V)  # data matrix
  C <- matrix(rnorm(n * T), n, T)  # trial regressors
  
  # run_fastlss
  res_a <- run_fastlss(A, Y, C)
  
  expect_true(is.matrix(res_a))
  expect_equal(dim(res_a), c(T, V))
})

test_that("run_stablss returns matrix with correct dims", {
  skip_if_not_installed("hrfals")
  
  # Use same test data
  n <- 100  
  q <- 3    
  V <- 4    
  
  A <- matrix(rnorm(n * q), n, q)
  Y <- matrix(rnorm(n * V), n, V)
  
  # run_stablss
  res_b <- run_stablss(A, Y)
  
  expect_true(is.matrix(res_b))
  expect_equal(dim(res_b), c(q, V))
})
