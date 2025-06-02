# Tests for Core Voxel-wise Fit Functions (Component 1)
# Tests for MHRF-CORE-VOXFIT-01 through MHRF-CORE-VOXFIT-05

test_that("project_out_confounds_core works without confounds", {
  # Create test data
  set.seed(123)
  n <- 50  # timepoints
  V <- 20  # voxels
  p <- 10  # HRF length
  k <- 3   # conditions
  
  Y_data <- matrix(rnorm(n * V), n, V)
  X_list <- lapply(1:k, function(i) matrix(rnorm(n * p), n, p))
  
  # Test with NULL confounds
  result <- project_out_confounds_core(Y_data, X_list, NULL)
  
  # Should return original matrices unchanged
  expect_identical(result$Y_proj_matrix, Y_data)
  expect_identical(result$X_list_proj_matrices, X_list)
})

test_that("project_out_confounds_core projects out confounds correctly", {
  # Create test data
  set.seed(456)
  n <- 100
  V <- 30
  p <- 15
  k <- 2
  q <- 4  # number of confounds
  
  # Create confounds (intercept + linear trend + quadratic + cubic)
  Z_confounds <- cbind(1, poly(1:n, degree = 3, raw = TRUE))
  
  # Create data with known confound contribution
  confound_effects_Y <- Z_confounds %*% matrix(rnorm(q * V), q, V)
  Y_clean <- matrix(rnorm(n * V), n, V)
  Y_data <- Y_clean + confound_effects_Y
  
  # Create design matrices with confound contribution
  X_list <- lapply(1:k, function(i) {
    X_clean <- matrix(rnorm(n * p), n, p)
    confound_effects_X <- Z_confounds %*% matrix(rnorm(q * p), q, p)
    X_clean + confound_effects_X
  })
  
  # Project out confounds
  result <- project_out_confounds_core(Y_data, X_list, Z_confounds)
  
  # Check dimensions
  expect_equal(dim(result$Y_proj_matrix), dim(Y_data))
  expect_equal(length(result$X_list_proj_matrices), k)
  for (i in 1:k) {
    expect_equal(dim(result$X_list_proj_matrices[[i]]), dim(X_list[[i]]))
  }
  
  # Check that projected data is orthogonal to confounds
  # The projection should remove all variance explained by confounds
  # Check that Z'Y_proj ≈ 0 relative to the scale of the data
  Y_proj_confound_corr <- crossprod(Z_confounds, result$Y_proj_matrix)
  Y_scale <- sqrt(sum(result$Y_proj_matrix^2))
  Z_scale <- sqrt(sum(Z_confounds^2))
  relative_error_Y <- sqrt(sum(Y_proj_confound_corr^2)) / (Y_scale * Z_scale / n)
  expect_lt(relative_error_Y, 1e-8)
  
  # Check that projected designs are orthogonal to confounds
  for (i in 1:k) {
    X_proj_confound_corr <- crossprod(Z_confounds, result$X_list_proj_matrices[[i]])
    X_scale <- sqrt(sum(result$X_list_proj_matrices[[i]]^2))
    relative_error_X <- sqrt(sum(X_proj_confound_corr^2)) / (X_scale * Z_scale / n)
    expect_lt(relative_error_X, 1e-7)
  }
})

test_that("project_out_confounds_core validates inputs correctly", {
  # Create valid base data
  n <- 50
  V <- 20
  p <- 10
  k <- 2
  
  Y_data <- matrix(rnorm(n * V), n, V)
  X_list <- lapply(1:k, function(i) matrix(rnorm(n * p), n, p))
  Z_confounds <- cbind(1, 1:n)
  
  # Test non-matrix Y
  expect_error(
    project_out_confounds_core(data.frame(Y_data), X_list, Z_confounds),
    "Y_data_matrix must be a matrix"
  )
  
  # Test non-list X
  expect_error(
    project_out_confounds_core(Y_data, X_list[[1]], Z_confounds),
    "Design matrices must be provided as a non-empty list"
  )
  
  # Test mismatched dimensions
  X_bad <- X_list
  X_bad[[1]] <- matrix(rnorm((n-1) * p), n-1, p)
  expect_error(
    project_out_confounds_core(Y_data, X_bad, Z_confounds),
    "Design matrix 1 has 49 rows but expected 50"
  )
  
  # Test too many confounds
  Z_bad <- matrix(rnorm(n * n), n, n)
  expect_error(
    project_out_confounds_core(Y_data, X_list, Z_bad),
    "too many columns"
  )
})

test_that("project_out_confounds_core handles rank deficient confounds with warning", {
  n <- 40
  V <- 10
  p <- 8
  k <- 2

  Y_data <- matrix(rnorm(n * V), n, V)
  X_list <- lapply(1:k, function(i) matrix(rnorm(n * p), n, p))

  # Create rank deficient confounds (second column duplicate)
  Z_confounds <- cbind(1:n, 2 * (1:n))

  expect_warning(
    res <- project_out_confounds_core(Y_data, X_list, Z_confounds),
    "rank deficient"
  )

  expect_equal(dim(res$Y_proj_matrix), dim(Y_data))
  expect_equal(length(res$X_list_proj_matrices), k)
})

test_that("project_out_confounds_core errors on NA confounds", {
  n <- 20
  V <- 5
  p <- 4
  k <- 1

  Y_data <- matrix(rnorm(n * V), n, V)
  X_list <- list(matrix(rnorm(n * p), n, p))
  Z_confounds <- cbind(1:n, rep(NA, n))

  expect_error(
    project_out_confounds_core(Y_data, X_list, Z_confounds),
    "non-finite values"
  )
})

test_that("project_out_confounds_core preserves data structure after projection", {
  # Create structured test data
  set.seed(789)
  n <- 80
  V <- 25
  p <- 12
  k <- 3
  
  # Create data with known signal not related to confounds
  t <- (1:n) / n * 2 * pi
  signal <- sin(2 * t) %*% t(runif(V))
  
  # Add confounds
  Z_confounds <- cbind(1, 1:n, (1:n)^2)
  confound_contrib <- Z_confounds %*% matrix(rnorm(3 * V), 3, V)
  Y_data <- signal + confound_contrib + 0.1 * matrix(rnorm(n * V), n, V)
  
  # Simple design matrices
  X_list <- lapply(1:k, function(i) {
    X <- matrix(0, n, p)
    # Put some structure in the design
    X[seq(i, n, by = k), ] <- matrix(rnorm(length(seq(i, n, by = k)) * p), ncol = p)
    X
  })
  
  # Project out confounds
  result <- project_out_confounds_core(Y_data, X_list, Z_confounds)
  
  # The signal structure should still be present (though not identical due to projection)
  # Check correlation between original signal and projected data
  signal_proj <- signal - Z_confounds %*% solve(crossprod(Z_confounds)) %*% crossprod(Z_confounds, signal)
  
  # Correlation should be high (signal preserved after projection)
  correlations <- diag(cor(signal_proj, result$Y_proj_matrix))
  expect_true(median(correlations) > 0.8)
})

test_that("transform_designs_to_manifold_basis_core works correctly", {
  # Create test data
  set.seed(123)
  n <- 100  # timepoints
  p <- 20   # HRF length
  m <- 5    # manifold dimensions
  k <- 3    # conditions
  
  # Create design matrices and reconstructor
  X_list <- lapply(1:k, function(i) matrix(rnorm(n * p), n, p))
  B_reconstructor <- matrix(rnorm(p * m), p, m)
  
  # Transform to manifold basis
  Z_list <- transform_designs_to_manifold_basis_core(X_list, B_reconstructor)
  
  # Check output structure
  expect_type(Z_list, "list")
  expect_length(Z_list, k)
  
  # Check dimensions of transformed matrices
  for (i in 1:k) {
    expect_equal(dim(Z_list[[i]]), c(n, m))
  }
  
  # Check that transformation is correct
  # Manually compute first transformation
  Z_expected_1 <- X_list[[1]] %*% B_reconstructor
  expect_equal(Z_list[[1]], Z_expected_1)
})

test_that("transform_designs_to_manifold_basis_core validates inputs", {
  # Valid inputs
  X_list <- list(matrix(1:20, 10, 2), matrix(21:40, 10, 2))
  B <- matrix(1:6, 2, 3)
  
  # Test non-list input
  expect_error(
    transform_designs_to_manifold_basis_core(matrix(1:20, 10, 2), B),
    "must be a list"
  )
  
  # Test non-matrix B
  expect_error(
    transform_designs_to_manifold_basis_core(X_list, data.frame(B)),
    "B_reconstructor_matrix must be a matrix"
  )
  
  # Test empty list
  expect_error(
    transform_designs_to_manifold_basis_core(list(), B),
    "cannot be empty"
  )
  
  # Test dimension mismatch
  X_bad <- list(matrix(1:30, 10, 3), matrix(1:20, 10, 2))  # Different p
  expect_error(
    transform_designs_to_manifold_basis_core(X_bad, B),
    "has 3 columns but B_reconstructor_matrix has 2 rows"
  )
})

test_that("transform_designs_to_manifold_basis_core preserves linear relationships", {
  # Test that linear combinations are preserved
  set.seed(456)
  n <- 50
  p <- 15
  m <- 4
  k <- 2
  
  # Create design matrices
  X_list <- lapply(1:k, function(i) matrix(rnorm(n * p), n, p))
  B_reconstructor <- matrix(rnorm(p * m), p, m)
  
  # Transform
  Z_list <- transform_designs_to_manifold_basis_core(X_list, B_reconstructor)
  
  # Linear combination of X's
  alpha <- c(0.3, 0.7)
  X_combined <- alpha[1] * X_list[[1]] + alpha[2] * X_list[[2]]
  
  # Should equal linear combination of Z's
  Z_combined_expected <- X_combined %*% B_reconstructor
  Z_combined_actual <- alpha[1] * Z_list[[1]] + alpha[2] * Z_list[[2]]
  
  expect_equal(Z_combined_actual, Z_combined_expected, tolerance = 1e-8)
})

test_that("solve_glm_for_gamma_core works correctly", {
  # Create test data
  set.seed(123)
  n <- 100  # timepoints
  m <- 4    # manifold dimensions
  k <- 3    # conditions
  V <- 50   # voxels
  
  # Create design matrices in manifold basis
  Z_list <- lapply(1:k, function(i) matrix(rnorm(n * m), n, m))
  
  # Create data with known structure
  # True gamma coefficients
  true_gamma <- matrix(rnorm((k * m) * V), k * m, V)
  X_tilde <- do.call(cbind, Z_list)
  Y_proj <- X_tilde %*% true_gamma + 0.1 * matrix(rnorm(n * V), n, V)
  
  # Solve GLM
  gamma_est <- solve_glm_for_gamma_core(Z_list, Y_proj, lambda_gamma = 0.01)
  
  # Check dimensions
  expect_equal(dim(gamma_est), c(k * m, V))
  
  # Check that estimates are close to truth (with some tolerance for noise and ridge)
  correlation <- cor(as.vector(true_gamma), as.vector(gamma_est))
  expect_gt(correlation, 0.9)
})

test_that("solve_glm_for_gamma_core orthogonal approximation works", {
  # Create test data with orthogonal conditions
  set.seed(456)
  n <- 150
  m <- 3
  k <- 4
  V <- 30
  
  # Create orthogonal design matrices
  Z_list <- list()
  for (i in 1:k) {
    Z <- matrix(0, n, m)
    # Make conditions temporally non-overlapping
    idx <- seq(from = (i-1) * (n/k) + 1, length.out = n/k)
    Z[idx, ] <- matrix(rnorm(length(idx) * m), length(idx), m)
    Z_list[[i]] <- Z
  }
  
  Y_proj <- matrix(rnorm(n * V), n, V)
  
  # Solve with and without orthogonal approximation
  gamma_full <- solve_glm_for_gamma_core(Z_list, Y_proj, lambda_gamma = 0.01, 
                                        orthogonal_approx_flag = FALSE)
  gamma_ortho <- solve_glm_for_gamma_core(Z_list, Y_proj, lambda_gamma = 0.01, 
                                          orthogonal_approx_flag = TRUE)
  
  # For truly orthogonal conditions, results should be very similar
  # But not identical due to ridge penalty structure
  correlation <- cor(as.vector(gamma_full), as.vector(gamma_ortho))
  expect_gt(correlation, 0.95)
  
  # Check that the approximation gives reasonable results
  expect_equal(dim(gamma_ortho), c(k * m, V))
})

test_that("solve_glm_for_gamma_core validates inputs", {
  # Valid base inputs
  Z_list <- list(matrix(1:20, 10, 2), matrix(21:40, 10, 2))
  Y <- matrix(1:50, 10, 5)
  
  # Test non-list Z
  expect_error(
    solve_glm_for_gamma_core(matrix(1:20, 10, 2), Y, 0.1),
    "must be a list"
  )
  
  # Test empty list
  expect_error(
    solve_glm_for_gamma_core(list(), Y, 0.1),
    "cannot be empty"
  )
  
  # Test non-matrix Y
  expect_error(
    solve_glm_for_gamma_core(Z_list, data.frame(Y), 0.1),
    "Y_proj_matrix must be a matrix"
  )
  
  # Test invalid lambda
  expect_error(
    solve_glm_for_gamma_core(Z_list, Y, -0.1),
    "non-negative scalar"
  )
  
  # Test dimension mismatch
  Z_bad <- list(matrix(1:20, 10, 2), matrix(1:18, 9, 2))
  expect_error(
    solve_glm_for_gamma_core(Z_bad, Y, 0.1),
    "must have 10 rows"
  )
  
  # Test inconsistent m dimensions
  Z_bad2 <- list(matrix(1:20, 10, 2), matrix(1:30, 10, 3))
  expect_error(
    solve_glm_for_gamma_core(Z_bad2, Y, 0.1),
    "same number of columns"
  )
})

test_that("solve_glm_for_gamma_core handles ridge regression correctly", {
  # Test that ridge penalty has expected effect
  set.seed(789)
  n <- 80
  m <- 5
  k <- 2
  V <- 20
  
  # Create slightly ill-conditioned problem
  Z1 <- matrix(rnorm(n * m), n, m)
  Z2 <- Z1 + 0.1 * matrix(rnorm(n * m), n, m)  # Nearly collinear
  Z_list <- list(Z1, Z2)
  
  Y_proj <- matrix(rnorm(n * V), n, V)
  
  # Solve with different ridge penalties
  gamma_small_ridge <- solve_glm_for_gamma_core(Z_list, Y_proj, lambda_gamma = 0.001)
  gamma_large_ridge <- solve_glm_for_gamma_core(Z_list, Y_proj, lambda_gamma = 1.0)
  
  # Larger ridge should give smaller coefficient magnitudes
  norm_small <- norm(gamma_small_ridge, "F")
  norm_large <- norm(gamma_large_ridge, "F")
  expect_lt(norm_large, norm_small)
  
  # Both should give valid results
  expect_false(any(is.na(gamma_small_ridge)))
  expect_false(any(is.na(gamma_large_ridge)))
})

test_that("extract_xi_beta_raw_svd_core works correctly", {
  # Create test data
  set.seed(123)
  m <- 4  # manifold dimensions
  k <- 3  # conditions
  V <- 50 # voxels
  
  # Create known Xi and Beta
  true_Xi <- matrix(rnorm(m * V), m, V)
  true_Beta <- matrix(rnorm(k * V), k, V)
  
  # Construct Gamma from Xi and Beta
  Gamma_coeffs <- matrix(0, k * m, V)
  for (vx in 1:V) {
    # Create k x m matrix: Beta (k×1) %*% Xi' (1×m)
    G_vx <- outer(true_Beta[, vx], true_Xi[, vx])
    # Vectorize as byrow=TRUE to match extract function
    Gamma_coeffs[, vx] <- as.vector(t(G_vx))
  }
  
  # Extract Xi and Beta using SVD
  result <- extract_xi_beta_raw_svd_core(Gamma_coeffs, m, k)
  
  # Check dimensions
  expect_equal(dim(result$Xi_raw_matrix), c(m, V))
  expect_equal(dim(result$Beta_raw_matrix), c(k, V))
  
  # Check reconstruction (up to sign)
  # Reconstruct Gamma from extracted Xi and Beta
  Gamma_reconstructed <- matrix(0, k * m, V)
  for (vx in 1:V) {
    # Beta (k×1) %*% Xi' (1×m) = k×m matrix
    G_vx_recon <- outer(result$Beta_raw_matrix[, vx], result$Xi_raw_matrix[, vx])
    # Vectorize as byrow to match input
    Gamma_reconstructed[, vx] <- as.vector(t(G_vx_recon))
  }
  
  # Should be very close (allowing for sign flips)
  reconstruction_error <- norm(abs(Gamma_coeffs) - abs(Gamma_reconstructed), "F") / norm(Gamma_coeffs, "F")
  expect_lt(reconstruction_error, 1e-5)  # Allow for numerical errors from robust SVD
})

test_that("extract_xi_beta_raw_svd_core handles near-zero singular values", {
  # Create test data with some zero columns
  set.seed(456)
  m <- 3
  k <- 2
  V <- 20
  
  # Create Gamma with some zero columns (voxels with no signal)
  Gamma_coeffs <- matrix(rnorm((k * m) * V), k * m, V)
  # Set some columns to zero
  zero_voxels <- c(5, 10, 15)
  Gamma_coeffs[, zero_voxels] <- 0
  
  # Extract Xi and Beta
  result <- extract_xi_beta_raw_svd_core(Gamma_coeffs, m, k)
  
  # Check that zero voxels have zero Xi and Beta
  expect_equal(result$Xi_raw_matrix[, zero_voxels], matrix(0, m, length(zero_voxels)))
  expect_equal(result$Beta_raw_matrix[, zero_voxels], matrix(0, k, length(zero_voxels)))
  
  # Non-zero voxels should have non-zero values
  non_zero_voxels <- setdiff(1:V, zero_voxels)
  expect_true(all(colSums(abs(result$Xi_raw_matrix[, non_zero_voxels])) > 0))
  expect_true(all(colSums(abs(result$Beta_raw_matrix[, non_zero_voxels])) > 0))
})

test_that("extract_xi_beta_raw_svd_core validates inputs", {
  # Valid gamma matrix
  gamma <- matrix(1:30, 6, 5)  # (2*3) x 5
  
  # Test non-matrix input
  expect_error(
    extract_xi_beta_raw_svd_core(data.frame(gamma), 2, 3),
    "must be a matrix"
  )
  
  # Test invalid m
  expect_error(
    extract_xi_beta_raw_svd_core(gamma, 0, 3),
    "positive integer"
  )
  
  expect_error(
    extract_xi_beta_raw_svd_core(gamma, 2.5, 3),
    "positive integer"
  )
  
  # Test invalid k
  expect_error(
    extract_xi_beta_raw_svd_core(gamma, 2, -1),
    "positive integer"
  )
  
  # Test dimension mismatch
  expect_error(
    extract_xi_beta_raw_svd_core(gamma, 3, 3),  # 3*3=9 != 6
    "has 6 rows but expected 9"
  )
})

test_that("extract_xi_beta_raw_svd_core preserves rank-1 structure", {
  # Test with exactly rank-1 data
  set.seed(789)
  m <- 5
  k <- 4
  V <- 30
  
  # Create rank-1 Gamma for each voxel
  Xi_rank1 <- matrix(rnorm(m * V), m, V)
  Beta_rank1 <- matrix(rnorm(k * V), k, V)
  
  Gamma_rank1 <- matrix(0, k * m, V)
  for (vx in 1:V) {
    # Create k x m matrix: Beta (k×1) %*% Xi' (1×m)
    G_vx <- Beta_rank1[, vx] %*% t(Xi_rank1[, vx])
    # Vectorize in row-major order (byrow=TRUE) to match extract function
    Gamma_rank1[, vx] <- as.vector(t(G_vx))
  }
  
  # Extract
  result <- extract_xi_beta_raw_svd_core(Gamma_rank1, m, k)
  
  # Check that we recover the same subspace (up to scaling and sign)
  for (vx in 1:V) {
    # The extract function reshapes gamma as k×m matrix
    gamma_v <- Gamma_rank1[, vx]
    G_reshaped <- matrix(gamma_v, k, m, byrow = TRUE)
    
    # The reconstruction from extracted components should match reshaped gamma
    # Beta (k×1) %*% Xi' (1×m) = k×m matrix
    G_extracted <- outer(result$Beta_raw_matrix[, vx], result$Xi_raw_matrix[, vx])
    
    # Should be identical up to numerical precision
    # Note: SVD reconstruction can introduce small numerical errors
    expect_equal(G_reshaped, G_extracted, tolerance = 1e-3)
  }
})

test_that("apply_intrinsic_identifiability_core works correctly", {
  # Create test data
  set.seed(123)
  m <- 4   # manifold dims
  k <- 3   # conditions
  V <- 50  # voxels
  p <- 20  # HRF length
  
  # Create test matrices
  Xi_raw <- matrix(rnorm(m * V), m, V)
  Beta_raw <- matrix(rnorm(k * V), k, V)
  B_reconstructor <- matrix(rnorm(p * m), p, m)
  
  # Create canonical HRF (simple gamma-like shape)
  t <- seq(0, p-1) / 2  # time in seconds
  h_ref <- dgamma(t, shape = 6, rate = 1) - 0.35 * dgamma(t, shape = 16, rate = 1)
  h_ref <- h_ref / max(h_ref)  # normalize to peak 1
  
  # Test with l2_norm scaling
  result <- apply_intrinsic_identifiability_core(
    Xi_raw, Beta_raw, B_reconstructor, h_ref,
    ident_scale_method = "l2_norm"
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("Xi_ident_matrix", "Beta_ident_matrix"))
  expect_equal(dim(result$Xi_ident_matrix), c(m, V))
  expect_equal(dim(result$Beta_ident_matrix), c(k, V))
  
  # Check that L2 norm constraint is satisfied
  for (vx in 1:V) {
    if (any(result$Xi_ident_matrix[, vx] != 0)) {
      hrf_vx <- B_reconstructor %*% result$Xi_ident_matrix[, vx]
      l2_norm <- sqrt(sum(hrf_vx^2))
      expect_equal(l2_norm, 1, tolerance = 1e-8)
    }
  }
})

test_that("apply_intrinsic_identifiability_core handles different scaling methods", {
  set.seed(456)
  m <- 3
  V <- 20
  p <- 15
  k <- 2
  
  Xi_raw <- matrix(rnorm(m * V), m, V)
  Beta_raw <- matrix(rnorm(k * V), k, V)
  B_reconstructor <- matrix(rnorm(p * m), p, m)
  h_ref <- rnorm(p)
  
  # Test all scaling methods
  result_l2 <- apply_intrinsic_identifiability_core(
    Xi_raw, Beta_raw, B_reconstructor, h_ref,
    ident_scale_method = "l2_norm"
  )
  
  result_max <- apply_intrinsic_identifiability_core(
    Xi_raw, Beta_raw, B_reconstructor, h_ref,
    ident_scale_method = "max_abs_val"
  )
  
  result_none <- apply_intrinsic_identifiability_core(
    Xi_raw, Beta_raw, B_reconstructor, h_ref,
    ident_scale_method = "none"
  )
  
  # Check scaling properties
  for (vx in 1:V) {
    if (any(Xi_raw[, vx] != 0)) {
      # L2 norm
      hrf_l2 <- B_reconstructor %*% result_l2$Xi_ident_matrix[, vx]
      expect_equal(sqrt(sum(hrf_l2^2)), 1, tolerance = 1e-8)
      
      # Max abs val
      hrf_max <- B_reconstructor %*% result_max$Xi_ident_matrix[, vx]
      expect_equal(max(abs(hrf_max)), 1, tolerance = 1e-8)
      
      # None - should preserve original scale times sign
      xi_orig <- Xi_raw[, vx]
      xi_none <- result_none$Xi_ident_matrix[, vx]
      # Check they differ only by sign
      expect_true(all(abs(abs(xi_orig) - abs(xi_none)) < 1e-10))
    }
  }
})

test_that("apply_intrinsic_identifiability_core preserves signal", {
  # Test that Xi * Beta product is preserved (up to sign)
  set.seed(789)
  m <- 5
  k <- 4
  V <- 30
  p <- 25
  
  Xi_raw <- matrix(rnorm(m * V), m, V)
  Beta_raw <- matrix(rnorm(k * V), k, V)
  B_reconstructor <- matrix(rnorm(p * m), p, m)
  h_ref <- rnorm(p)
  
  # Apply identifiability
  result <- apply_intrinsic_identifiability_core(
    Xi_raw, Beta_raw, B_reconstructor, h_ref
  )
  
  # Check signal preservation for each voxel
  for (vx in 1:V) {
    for (cond in 1:k) {
      # Original signal contribution
      signal_orig <- sum(Xi_raw[, vx] * Beta_raw[cond, vx])
      
      # Identifiable signal contribution
      signal_ident <- sum(result$Xi_ident_matrix[, vx] * result$Beta_ident_matrix[cond, vx])
      
      # Should be equal up to sign
      expect_equal(abs(signal_orig), abs(signal_ident), tolerance = 1e-8)
    }
  }
})

test_that("apply_intrinsic_identifiability_core validates inputs", {
  # Valid inputs
  Xi <- matrix(1:10, 2, 5)
  Beta <- matrix(1:15, 3, 5)
  B <- matrix(1:8, 4, 2)
  h_ref <- 1:4
  
  # Test matrix inputs
  expect_error(
    apply_intrinsic_identifiability_core(data.frame(Xi), Beta, B, h_ref),
    "Xi_raw_matrix must be a matrix"
  )
  
  expect_error(
    apply_intrinsic_identifiability_core(Xi, data.frame(Beta), B, h_ref),
    "Beta_raw_matrix must be a matrix"
  )
  
  expect_error(
    apply_intrinsic_identifiability_core(Xi, Beta, data.frame(B), h_ref),
    "B_reconstructor_matrix must be a matrix"
  )
  
  # Test h_ref vector
  expect_error(
    apply_intrinsic_identifiability_core(Xi, Beta, B, matrix(h_ref)),
    "h_ref_shape_vector must be a numeric vector"
  )
  
  # Test dimension mismatches
  expect_error(
    apply_intrinsic_identifiability_core(Xi, matrix(1:12, 3, 4), B, h_ref),
    "must have the same number of columns"
  )
  
  expect_error(
    apply_intrinsic_identifiability_core(Xi, Beta, matrix(1:12, 4, 3), h_ref),
    "must have m columns"
  )
  
  expect_error(
    apply_intrinsic_identifiability_core(Xi, Beta, B, 1:5),
    "must have length p"
  )
  
  # Test invalid methods
  expect_error(
    apply_intrinsic_identifiability_core(Xi, Beta, B, h_ref, ident_scale_method = "invalid"),
    "ident_scale_method must be one of"
  )
  
  expect_error(
    apply_intrinsic_identifiability_core(Xi, Beta, B, h_ref, ident_sign_method = "invalid"),
    "ident_sign_method must be one of"
  )
})

test_that("apply_intrinsic_identifiability_core handles zero voxels", {
  # Test with some zero Xi columns
  set.seed(321)
  m <- 3
  k <- 2
  V <- 10
  p <- 12
  
  Xi_raw <- matrix(rnorm(m * V), m, V)
  Beta_raw <- matrix(rnorm(k * V), k, V)
  
  # Set some voxels to zero
  zero_voxels <- c(2, 5, 8)
  Xi_raw[, zero_voxels] <- 0
  
  B_reconstructor <- matrix(rnorm(p * m), p, m)
  h_ref <- rnorm(p)
  
  # Apply identifiability
  result <- apply_intrinsic_identifiability_core(
    Xi_raw, Beta_raw, B_reconstructor, h_ref
  )
  
  # Zero voxels should remain zero
  expect_equal(result$Xi_ident_matrix[, zero_voxels], matrix(0, m, length(zero_voxels)))
  
  # Both Xi and Beta should be zero for these voxels
  expect_equal(result$Beta_ident_matrix[, zero_voxels], matrix(0, k, length(zero_voxels)))
})

test_that("apply_intrinsic_identifiability_core zeros tiny HRFs", {
  set.seed(987)
  m <- 3
  k <- 2
  V <- 5
  p <- 10

  Xi_raw <- matrix(rnorm(m * V), m, V)
  Beta_raw <- matrix(rnorm(k * V), k, V)
  B_reconstructor <- matrix(rnorm(p * m), p, m)
  h_ref <- rnorm(p)

  tiny_voxel <- 1
  Xi_raw[, tiny_voxel] <- 1e-12

  result <- apply_intrinsic_identifiability_core(
    Xi_raw, Beta_raw, B_reconstructor, h_ref,
    ident_scale_method = "l2_norm"
  )

  expect_equal(result$Xi_ident_matrix[, tiny_voxel], rep(0, m))
  expect_equal(result$Beta_ident_matrix[, tiny_voxel], rep(0, k))
})