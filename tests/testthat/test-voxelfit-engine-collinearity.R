context("core voxel fit engine stability")
library(testthat)
library(manifoldhrf)

# Perfect Collinearity

test_that("solve_glm_for_gamma_core errors on perfectly collinear design", {
  set.seed(1)
  n <- 20
  m <- 1
  V <- 1
  # Two conditions with identical design columns
  base <- rnorm(n)
  Z1 <- matrix(base, n, m)
  Z2 <- Z1
  Y <- matrix(rnorm(n * V), n, V)
  expect_error(
    solve_glm_for_gamma_core(list(Z1, Z2), Y, lambda_gamma = 0),
    "singular"
  )
})

# Near Collinearity

test_that("solve_glm_for_gamma_core produces unstable betas for near collinear design", {
  set.seed(2)
  n <- 30
  m <- 1
  V <- 1
  base <- rnorm(n)
  Z1 <- matrix(base, n, m)
  Z2 <- matrix(base + rnorm(n, sd = 1e-4), n, m)
  Y <- matrix(rnorm(n * V), n, V)
  gamma_est <- solve_glm_for_gamma_core(list(Z1, Z2), Y, lambda_gamma = 0)
  expect_equal(dim(gamma_est), c(2 * m, V))
  expect_true(max(abs(gamma_est)) > 100)
})

# Zero Variance Voxels

test_that("estimate_final_condition_betas_core handles zero variance voxels", {
  set.seed(3)
  n <- 10
  V <- 5
  p <- 3
  k <- 2
  Y <- matrix(0, n, V)
  Xc <- list(matrix(rnorm(n * p), n, p), matrix(rnorm(n * p), n, p))
  H <- matrix(rnorm(p * V), p, V)
  res <- estimate_final_condition_betas_core(Y, Xc, H)
  expect_equal(dim(res), c(k, V))
  expect_true(all(is.finite(res)))
  expect_true(all(res == 0))
})

