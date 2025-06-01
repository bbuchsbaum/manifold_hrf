# Tests for mhrf_lss_parameters helper

test_that("mhrf_lss_parameters merges preset and user values", {
  base <- mhrf_lss_parameters("balanced")
  override <- mhrf_lss_parameters("balanced", lambda_gamma = 0.05, use_parallel = TRUE)

  expect_true(is.list(base))
  expect_equal(override$lambda_gamma, 0.05)
  expect_true(override$use_parallel)
  expect_equal(base$m_manifold_dim_target, override$m_manifold_dim_target)
})
