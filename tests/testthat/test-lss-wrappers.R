library(testthat)

skip_if_not_installed("hrfals")

# simple matrices for testing
n <- 100  # timepoints
q <- 3    # nuisance regressors
V <- 4    # voxels  
T <- 5    # trials

A <- matrix(rnorm(n * q), n, q)  # design matrix
Y <- matrix(rnorm(n * V), n, V)  # data matrix
C <- matrix(rnorm(n * T), n, T)  # trial regressors

# expect that hrfals functions return matrices of dim T x V

# run_fastlss
res_a <- run_fastlss(A, Y, C)

test_that("run_fastlss returns matrix with correct dims", {
  expect_true(is.matrix(res_a))
  expect_equal(dim(res_a), c(T, V))
})

# run_stablss
res_b <- run_stablss(A, Y)

test_that("run_stablss returns matrix with correct dims", {
  expect_true(is.matrix(res_b))
  expect_equal(dim(res_b), c(q, V))
})
