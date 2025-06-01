context("fallback cascade")

 test_that("run_with_fallback_cascade falls back when manifold fails", {
  set.seed(1)
  Y <- matrix(rnorm(20), 5, 4)
  X <- list(matrix(rnorm(5), 5, 1))
  params <- list(TR = 1, m_manifold_dim = -2)
  result <- run_with_fallback_cascade(Y, X, params)
  expect_true(result$success)
  expect_true(all(result$method == "pca"))
  expect_true(is.matrix(result$results$Beta))
 })
