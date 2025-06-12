context("core voxel fit engine")
library(testthat)
library(manifoldhrf)

set.seed(1)

m <- 2
k <- 3
V <- 4
Gamma <- matrix(rnorm(m * k * V), m * k, V)

 test_that("extract_xi_beta_raw_svd_core returns correct dimensions", {
   res <- extract_xi_beta_raw_svd_core(Gamma, m, k)
   expect_equal(dim(res$Xi_raw_matrix), c(m, V))
   expect_equal(dim(res$Beta_raw_matrix), c(k, V))
 })

test_that("extract_xi_beta_raw_svd_core handles matrix orientation correctly", {
  # Test with specific values to verify correct reshape
  m_test <- 3  # manifold dimensions
  k_test <- 2  # conditions
  V_test <- 1  # single voxel for clarity
  
  # Create gamma vector ordered as [cond1_dim1, cond1_dim2, cond1_dim3, cond2_dim1, cond2_dim2, cond2_dim3]
  # This represents condition 1 = [1, 2, 3] and condition 2 = [4, 5, 6]
  gamma_vec <- 1:6
  Gamma_test <- matrix(gamma_vec, nrow = k_test * m_test, ncol = V_test)
  
  # Expected k x m matrix after reshape:
  # Row 1: condition 1 = [1, 2, 3]
  # Row 2: condition 2 = [4, 5, 6]
  expected_G <- matrix(c(1, 2, 3, 4, 5, 6), nrow = k_test, ncol = m_test, byrow = TRUE)
  
  # Compute SVD of expected matrix
  expected_svd <- svd(expected_G)
  
  # Run extraction
  res <- extract_xi_beta_raw_svd_core(Gamma_test, m_test, k_test)
  
  # Verify dimensions
  expect_equal(dim(res$Xi_raw_matrix), c(m_test, V_test))
  expect_equal(dim(res$Beta_raw_matrix), c(k_test, V_test))
  
  # Verify the SVD decomposition is correct
  # Reconstruct the matrix from extracted components
  xi_extracted <- res$Xi_raw_matrix[, 1]
  beta_extracted <- res$Beta_raw_matrix[, 1]
  
  # The rank-1 approximation should match the first singular value decomposition
  reconstructed <- outer(beta_extracted, xi_extracted)
  expected_rank1 <- expected_svd$d[1] * outer(expected_svd$u[, 1], expected_svd$v[, 1])
  
  expect_equal(reconstructed, expected_rank1, tolerance = 1e-10)
})

 test_that("apply_intrinsic_identifiability_core works", {
   Xi_raw <- matrix(rnorm(m * V), m, V)
   Beta_raw <- matrix(rnorm(k * V), k, V)
   B <- matrix(rnorm(5 * m), 5, m)
   h_ref <- rnorm(5)
   res <- apply_intrinsic_identifiability_core(Xi_raw, Beta_raw, B, h_ref)
   expect_equal(dim(res$Xi_ident_matrix), c(m, V))
   expect_equal(dim(res$Beta_ident_matrix), c(k, V))
 })

 test_that("make_voxel_graph_laplacian_core returns sparse Laplacian", {
   coords <- matrix(seq_len(9), ncol = 3)
   L <- make_voxel_graph_laplacian_core(coords, num_neighbors_Lsp = 2)
   expect_s4_class(L, "dgCMatrix")
   expect_equal(dim(L), c(nrow(coords), nrow(coords)))
 })

test_that("make_voxel_graph_laplacian_core produces correct Laplacian for 2-voxel chain", {
  coords <- matrix(c(0, 0, 0,
                     1, 0, 0), nrow = 2, byrow = TRUE)
  L <- make_voxel_graph_laplacian_core(coords, num_neighbors_Lsp = 1)
  Lm <- as.matrix(L)
  expect_equal(Lm, matrix(c(1, -1,
                             -1, 1), nrow = 2, byrow = TRUE))
})

test_that("apply_spatial_smoothing_core returns same dimension", {
  coords <- matrix(seq_len(9), ncol = 3)
  L <- make_voxel_graph_laplacian_core(coords, num_neighbors_Lsp = 2)
  Xi <- matrix(rnorm(m * nrow(coords)), m, nrow(coords))
  Xi_s <- apply_spatial_smoothing_core(Xi, L, 0.1)
  expect_equal(dim(Xi_s), dim(Xi))
})

test_that("prepare_lss_fixed_components_core returns matrices of correct size", {
  # Create a matrix where q_lss < n (3 columns, 10 rows)
  A <- cbind(1, matrix(rnorm(20), nrow = 10, ncol = 2))
  res <- prepare_lss_fixed_components_core(
    A_fixed_regressors_matrix = A, 
    lambda_ridge_A = 0.01
  )
  expect_null(res$P_lss)  # fmrilss handles internally
  expect_true(res$has_intercept)  # Should detect intercept column
})

test_that("reconstruct_hrf_shapes_core multiplies matrices", {
  B <- matrix(rnorm(10), 5, 2)
  Xi <- matrix(rnorm(2 * 3), 2, 3)
  H <- reconstruct_hrf_shapes_core(B, Xi)
  expect_equal(dim(H), c(5, 3))
})

test_that("run_lss_for_voxel returns vector of length T", {
  Y <- rnorm(5)
  X_list <- list(matrix(1:5, ncol = 1), matrix(5:1, ncol = 1))
  H <- rnorm(1)
  res <- run_lss_for_voxel(
    y_voxel = Y,
    X_trial_list = X_list,
    h_voxel = H,
    TR = 2
  )
  expect_length(res, length(X_list))
})

test_that("estimate_final_condition_betas_core returns matrix of correct dims", {
  Y <- matrix(rnorm(15), 5, 3)
  Xc <- list(matrix(1:5, ncol = 1), matrix(5:1, ncol = 1))
  H <- matrix(rnorm(1 * 3), 1, 3)
  res <- estimate_final_condition_betas_core(Y, Xc, H)
  expect_equal(dim(res), c(length(Xc), ncol(Y)))
})

