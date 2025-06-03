context("identifiability fallback")

set.seed(42)

# Simple identity reconstructor for clarity
B <- diag(3)
# reference HRF with mean 0 for zero correlation test
h_ref <- c(1, -1, 0)

# Xi that will produce HRF with zero correlation to h_ref and negative sum
# Need sum != 0 for sign to work properly
Xi_raw <- matrix(c(-1.5, -1.5, 2), nrow = 3, ncol = 1)
Beta_raw <- matrix(1, nrow = 1, ncol = 1)

# provide projected data and designs even though they should be ignored
Y_proj <- matrix(rnorm(3), 3, 1)
X_list <- list(matrix(rnorm(3 * 3), 3, 3))

res <- apply_intrinsic_identifiability_core(
  Xi_raw_matrix = Xi_raw,
  Beta_raw_matrix = Beta_raw,
  B_reconstructor_matrix = B,
  h_ref_shape_vector = h_ref,
  Y_proj_matrix = Y_proj,
  X_condition_list_proj_matrices = X_list,
  ident_sign_method = "canonical_correlation"
)

# Correlation with h_ref is zero so RMS rule should flip sign

test_that("canonical correlation falls back to RMS rule", {
  expect_equal(res$Xi_ident_matrix[1, 1] > 0, TRUE)
  expect_equal(res$Beta_ident_matrix[1, 1] < 0, TRUE)
})

# Second test: Xi that produces positive sum HRF with zero correlation
Xi_raw2 <- matrix(c(1, 1, -2), nrow = 3, ncol = 1)
Beta_raw2 <- matrix(1, nrow = 1, ncol = 1)
res2 <- apply_intrinsic_identifiability_core(
  Xi_raw_matrix = Xi_raw2,
  Beta_raw_matrix = Beta_raw2,
  B_reconstructor_matrix = B,
  h_ref_shape_vector = h_ref,
  Y_proj_matrix = Y_proj,
  X_condition_list_proj_matrices = X_list,
  ident_sign_method = "canonical_correlation"
)

test_that("RMS fallback preserves positive sign", {
  expect_equal(res2$Xi_ident_matrix[1, 1] > 0, TRUE)
  expect_equal(res2$Beta_ident_matrix[1, 1] > 0, TRUE)
})
