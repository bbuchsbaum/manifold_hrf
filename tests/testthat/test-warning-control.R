test_that("verbose_warnings parameter exists and function works", {
  # Create a simple test case
  set.seed(42)
  k <- 2  # conditions
  m <- 2  # manifold dim  
  V <- 3  # voxels
  
  # Create normal gamma matrix
  Gamma <- matrix(rnorm(k * m * V), nrow = k * m, ncol = V)
  
  # Test that function works with verbose_warnings = FALSE (default)
  result_quiet <- extract_xi_beta_raw_svd_core(
    Gamma_coeffs_matrix = Gamma,
    m_manifold_dim = m,
    k_conditions = k,
    method = "robust",
    verbose_warnings = FALSE
  )
  
  # Test that function works with verbose_warnings = TRUE
  result_verbose <- extract_xi_beta_raw_svd_core(
    Gamma_coeffs_matrix = Gamma,
    m_manifold_dim = m,
    k_conditions = k,
    method = "robust",
    verbose_warnings = TRUE
  )
  
  # Results should be identical regardless of warning setting
  expect_equal(result_quiet$Xi_raw_matrix, result_verbose$Xi_raw_matrix)
  expect_equal(result_quiet$Beta_raw_matrix, result_verbose$Beta_raw_matrix)
  
  # Check that both results have correct dimensions
  expect_equal(dim(result_quiet$Xi_raw_matrix), c(m, V))
  expect_equal(dim(result_quiet$Beta_raw_matrix), c(k, V))
  expect_equal(dim(result_verbose$Xi_raw_matrix), c(m, V))
  expect_equal(dim(result_verbose$Beta_raw_matrix), c(k, V))
}) 