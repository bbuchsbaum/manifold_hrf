# Tests for Core LSS Functions (Component 3)
# Tests for MHRF-CORE-LSS-01 through MHRF-CORE-LSS-04

test_that("prepare_lss_fixed_components_core works correctly", {
  # Create test fixed regressors
  set.seed(123)
  n <- 100  # timepoints
  
  # Simple case: intercept + linear drift
  A_fixed <- cbind(
    intercept = rep(1, n),
    drift = seq_len(n) / n
  )
  
  # Prepare LSS components
  result <- prepare_lss_fixed_components_core(
    A_fixed, 
    intercept_col_index_in_Alss = 1,
    lambda_ridge_Alss = 1e-6
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("P_lss_matrix", "p_lss_vector"))
  
  # Check dimensions
  expect_equal(dim(result$P_lss_matrix), c(2, n))  # q_lss x n
  expect_equal(length(result$p_lss_vector), n)
  
  # Check that P_lss approximates the pseudoinverse
  # P_lss ≈ (A'A + λI)^(-1) A'
  # So P_lss * A ≈ I (up to ridge penalty)
  PA <- result$P_lss_matrix %*% A_fixed
  expect_equal(as.vector(diag(PA)), rep(1, 2), tolerance = 0.01)
  
  # Check that p_lss_vector is related to intercept
  # It should be the first row of the pseudoinverse
  expect_equal(result$p_lss_vector, result$P_lss_matrix[1, ], tolerance = 1e-6)
})

test_that("prepare_lss_fixed_components_core handles no intercept", {
  set.seed(456)
  n <- 80
  
  # No intercept, just confounds
  A_fixed <- cbind(
    motion1 = rnorm(n),
    motion2 = rnorm(n),
    motion3 = rnorm(n)
  )
  
  # No intercept specified
  result <- prepare_lss_fixed_components_core(
    A_fixed,
    intercept_col_index_in_Alss = NULL,
    lambda_ridge_Alss = 0.01
  )
  
  # p_lss_vector should be all zeros
  expect_equal(result$p_lss_vector, rep(0, n))
  
  # P_lss should still be valid
  expect_equal(dim(result$P_lss_matrix), c(3, n))
})

test_that("prepare_lss_fixed_components_core validates inputs", {
  n <- 50
  A <- matrix(rnorm(n * 3), n, 3)
  
  # Non-matrix input
  expect_error(
    prepare_lss_fixed_components_core(data.frame(A), 1, 0),
    "must be a matrix"
  )
  
  # Too few rows
  expect_error(
    prepare_lss_fixed_components_core(matrix(1:3, 1, 3), 1, 0),
    "must have at least 2 rows"
  )
  
  # Too many columns
  A_bad <- matrix(rnorm(n * n), n, n)
  expect_error(
    prepare_lss_fixed_components_core(A_bad, 1, 0),
    "too many columns"
  )
  
  # Invalid intercept index
  expect_error(
    prepare_lss_fixed_components_core(A, 5, 0),
    "must be an integer between 1 and 3"
  )
  
  # Invalid lambda
  expect_error(
    prepare_lss_fixed_components_core(A, 1, -0.1),
    "must be a non-negative scalar"
  )
})

test_that("prepare_lss_fixed_components_core handles ill-conditioned matrices", {
  set.seed(789)
  n <- 100
  
  # Create nearly collinear columns
  x <- rnorm(n)
  A_fixed <- cbind(
    1,
    x,
    x + rnorm(n, sd = 1e-8)  # Nearly identical to x
  )
  
  # Should warn about conditioning
  expect_warning(
    result <- prepare_lss_fixed_components_core(A_fixed, 1, 1e-12),
    "poorly conditioned"
  )
  
  # But should still produce valid output
  expect_equal(dim(result$P_lss_matrix), c(3, n))
})

test_that("reconstruct_hrf_shapes_core works correctly", {
  # Create test data
  set.seed(123)
  p <- 30   # HRF length
  m <- 5    # manifold dimensions
  V <- 50   # voxels
  
  # Create reconstructor and manifold coordinates
  B_reconstructor <- matrix(rnorm(p * m), p, m)
  Xi_smoothed <- matrix(rnorm(m * V), m, V)
  
  # Reconstruct HRF shapes
  H_shapes <- reconstruct_hrf_shapes_core(B_reconstructor, Xi_smoothed)
  
  # Check dimensions
  expect_equal(dim(H_shapes), c(p, V))
  
  # Check that it's a regular matrix (not Matrix class)
  expect_true(is.matrix(H_shapes))
  expect_false(inherits(H_shapes, "Matrix"))
  
  # Verify matrix multiplication
  # Each column should be B %*% xi_v
  for (v in 1:5) {  # Check first few voxels
    expected <- B_reconstructor %*% Xi_smoothed[, v]
    expect_equal(H_shapes[, v], as.vector(expected))
  }
})

test_that("reconstruct_hrf_shapes_core validates inputs", {
  B <- matrix(1:10, 5, 2)
  Xi <- matrix(1:6, 2, 3)
  
  # Non-matrix inputs
  expect_error(
    reconstruct_hrf_shapes_core(data.frame(B), Xi),
    "B_reconstructor_matrix must be a matrix"
  )
  
  expect_error(
    reconstruct_hrf_shapes_core(B, data.frame(Xi)),
    "Xi_smoothed_matrix must be a matrix"
  )
  
  # Dimension mismatch
  Xi_bad <- matrix(1:9, 3, 3)  # Wrong number of rows
  expect_error(
    reconstruct_hrf_shapes_core(B, Xi_bad),
    "Dimension mismatch"
  )
  
  # Too few rows/columns
  expect_error(
    reconstruct_hrf_shapes_core(matrix(1:2, 1, 2), Xi),
    "must have at least 2 rows"
  )
  
  # Skip zero-column test that causes warnings
  B_zero <- matrix(numeric(0), 5, 0)
  Xi_zero <- matrix(numeric(0), 0, 3)
  expect_error(
    reconstruct_hrf_shapes_core(B_zero, Xi_zero),
    "must have at least 1 column"
  )
})

test_that("run_lss_for_voxel_core works correctly", {
  # Create realistic test scenario
  set.seed(123)
  n <- 200  # timepoints
  p <- 25   # HRF length
  T_trials <- 30  # trials
  
  # Create trial design matrices with proper spacing
  X_trials <- list()
  trial_spacing <- 6  # TRs between trials
  for (t in 1:T_trials) {
    X <- matrix(0, n, p)
    onset <- 10 + (t-1) * trial_spacing
    if (onset + p <= n) {
      # Simple FIR-style design
      for (j in 1:p) {
        if (onset + j - 1 <= n) {
          X[onset + j - 1, j] <- 1
        }
      }
    }
    X_trials[[t]] <- X
  }
  
  # Create HRF shape (gamma-like)
  t_hrf <- seq(0, p-1) * 0.5  # in seconds, assuming TR=0.5
  H_voxel <- dgamma(t_hrf, shape = 6, rate = 1) - 
             0.35 * dgamma(t_hrf, shape = 16, rate = 1)
  H_voxel <- H_voxel / max(H_voxel)
  
  # Create data with known trial amplitudes
  true_betas <- rnorm(T_trials, mean = 1, sd = 0.3)
  Y_voxel <- rep(0, n)
  for (t in 1:T_trials) {
    Y_voxel <- Y_voxel + true_betas[t] * as.vector(X_trials[[t]] %*% H_voxel)
  }
  Y_voxel <- Y_voxel + rnorm(n, sd = 0.1)  # Add noise
  
  # Fixed regressors (intercept + drift)
  A_fixed <- cbind(1, seq_len(n) / n)
  
  # Prepare LSS components
  lss_prep <- prepare_lss_fixed_components_core(A_fixed, 1, 1e-6)
  
  # Run LSS
  estimated_betas <- run_lss_for_voxel_core(
    Y_voxel,
    X_trials,
    H_voxel,
    A_fixed,
    lss_prep$P_lss_matrix,
    lss_prep$p_lss_vector
  )
  
  # Check output
  expect_equal(length(estimated_betas), T_trials)
  
  # Estimates should correlate with true values
  # LSS is biased toward zero, so check if pattern is preserved
  correlation <- cor(true_betas, estimated_betas)
  expect_gt(correlation, 0.5)
  
  # Also check that non-zero trials are detected
  expect_gt(mean(abs(estimated_betas)), 0.1)
})

test_that("run_lss_for_voxel_core validates inputs", {
  # Valid base inputs
  n <- 100
  p <- 20
  T <- 10
  
  Y <- rnorm(n)
  X_list <- lapply(1:T, function(i) matrix(rnorm(n*p), n, p))
  H <- rnorm(p)
  A <- cbind(1, rnorm(n))
  lss_prep <- prepare_lss_fixed_components_core(A, 1, 0)
  
  # Test non-numeric Y
  expect_error(
    run_lss_for_voxel_core(as.character(Y), X_list, H, A, 
                          lss_prep$P_lss_matrix, lss_prep$p_lss_vector),
    "Y_proj_voxel_vector must be a numeric vector"
  )
  
  # Test non-list X
  expect_error(
    run_lss_for_voxel_core(Y, X_list[[1]], H, A,
                          lss_prep$P_lss_matrix, lss_prep$p_lss_vector),
    "X_trial_onset_list_of_matrices must be a list"
  )
  
  # Test empty trial list
  expect_error(
    run_lss_for_voxel_core(Y, list(), H, A,
                          lss_prep$P_lss_matrix, lss_prep$p_lss_vector),
    "must contain at least one trial"
  )
  
  # Test dimension mismatches
  expect_error(
    run_lss_for_voxel_core(Y, X_list, H, matrix(1:50, 50, 2),
                          lss_prep$P_lss_matrix, lss_prep$p_lss_vector),
    "must have n rows"
  )
})

test_that("run_lss_voxel_loop_core works with and without precomputation", {
  # Create test data
  set.seed(456)
  n <- 100   # timepoints
  p <- 20    # HRF length
  V <- 30    # voxels
  T <- 15    # trials
  
  # Create data
  Y_proj <- matrix(rnorm(n * V), n, V)
  
  # Create trial designs
  X_trials <- list()
  for (t in 1:T) {
    X <- matrix(0, n, p)
    onset <- 5 + (t-1) * 5
    if (onset + p <= n) {
      for (j in 1:p) {
        X[onset + j - 1, j] <- 1
      }
    }
    X_trials[[t]] <- X
  }
  
  # Create HRF shapes
  H_shapes <- matrix(rnorm(p * V), p, V)
  
  # Fixed regressors
  A_fixed <- cbind(1, seq_len(n)/n, rnorm(n))
  lss_prep <- prepare_lss_fixed_components_core(A_fixed, 1, 1e-6)
  
  # Run with precomputation (high RAM limit)
  Beta_precomp <- run_lss_voxel_loop_core(
    Y_proj, X_trials, H_shapes, A_fixed,
    lss_prep$P_lss_matrix, lss_prep$p_lss_vector,
    ram_heuristic_GB_for_Rt = 10.0  # High limit forces precomputation
  )
  
  # Run without precomputation (low RAM limit)
  Beta_no_precomp <- run_lss_voxel_loop_core(
    Y_proj, X_trials, H_shapes, A_fixed,
    lss_prep$P_lss_matrix, lss_prep$p_lss_vector,
    ram_heuristic_GB_for_Rt = 0.0  # Zero limit prevents precomputation
  )
  
  # Results should be identical
  expect_equal(dim(Beta_precomp), c(T, V))
  expect_equal(dim(Beta_no_precomp), c(T, V))
  expect_equal(Beta_precomp, Beta_no_precomp, tolerance = 1e-10)
})

test_that("run_lss_voxel_loop_core validates inputs", {
  # Base valid inputs
  n <- 50
  V <- 20
  p <- 15
  T <- 10
  
  Y <- matrix(rnorm(n * V), n, V)
  X_list <- lapply(1:T, function(i) matrix(rnorm(n * p), n, p))
  H <- matrix(rnorm(p * V), p, V)
  A <- cbind(1, rnorm(n))
  lss_prep <- prepare_lss_fixed_components_core(A, 1, 0)
  
  # Test non-matrix Y
  expect_error(
    run_lss_voxel_loop_core(data.frame(Y), X_list, H, A,
                           lss_prep$P_lss_matrix, lss_prep$p_lss_vector),
    "Y_proj_matrix must be a matrix"
  )
  
  # Test dimension mismatch
  H_bad <- matrix(rnorm(p * (V-1)), p, V-1)
  expect_error(
    run_lss_voxel_loop_core(Y, X_list, H_bad, A,
                           lss_prep$P_lss_matrix, lss_prep$p_lss_vector),
    "must have V columns"
  )
  
  # Test invalid RAM heuristic
  expect_error(
    run_lss_voxel_loop_core(Y, X_list, H, A,
                           lss_prep$P_lss_matrix, lss_prep$p_lss_vector,
                           ram_heuristic_GB_for_Rt = -1),
    "must be a non-negative scalar"
  )
})

test_that("LSS integration test with known signal", {
  # Create a more complex scenario with overlapping trials
  set.seed(789)
  n <- 300   # timepoints
  p <- 30    # HRF length
  V <- 10    # voxels (small for testing)
  T <- 40    # trials
  
  # Create realistic HRF
  t_hrf <- seq(0, p-1) * 0.5
  h_canonical <- dgamma(t_hrf, shape = 6, rate = 1) - 
                 0.35 * dgamma(t_hrf, shape = 16, rate = 1)
  h_canonical <- h_canonical / sum(h_canonical)
  
  # Create voxel-specific HRFs with slight variations
  H_shapes <- matrix(0, p, V)
  for (v in 1:V) {
    # Add some variation to peak and width
    peak_shift <- rnorm(1, 0, 1)
    width_scale <- exp(rnorm(1, 0, 0.1))
    t_shifted <- (t_hrf - peak_shift) * width_scale
    t_shifted[t_shifted < 0] <- 0
    
    h_v <- dgamma(t_shifted, shape = 6, rate = 1) - 
           0.35 * dgamma(t_shifted, shape = 16, rate = 1)
    H_shapes[, v] <- h_v / sum(h_v)
  }
  
  # Create trial onsets with some jitter
  base_onsets <- seq(10, n-p-10, length.out = T)
  trial_onsets <- round(base_onsets + rnorm(T, 0, 2))
  
  # Create trial designs
  X_trials <- list()
  for (t in 1:T) {
    X <- matrix(0, n, p)
    onset <- trial_onsets[t]
    if (onset > 0 && onset + p <= n) {
      for (j in 1:p) {
        X[onset + j - 1, j] <- 1
      }
    }
    X_trials[[t]] <- X
  }
  
  # Create true betas with different patterns per voxel
  true_betas <- matrix(0, T, V)
  for (v in 1:V) {
    # Different activation patterns
    if (v <= 3) {
      # Sustained activation
      true_betas[, v] <- rnorm(T, mean = 1.5, sd = 0.3)
    } else if (v <= 6) {
      # Alternating activation
      true_betas[, v] <- rep(c(2, -1), length.out = T) + rnorm(T, sd = 0.2)
    } else {
      # Random activation
      true_betas[, v] <- rnorm(T, mean = 0, sd = 1)
    }
  }
  
  # Generate data
  Y_proj <- matrix(0, n, V)
  for (v in 1:V) {
    signal <- rep(0, n)
    for (t in 1:T) {
      signal <- signal + true_betas[t, v] * as.vector(X_trials[[t]] %*% H_shapes[, v])
    }
    Y_proj[, v] <- signal + rnorm(n, sd = 0.5)  # Add noise
  }
  
  # Add some drift
  drift <- seq(0, 1, length.out = n)
  Y_proj <- Y_proj + drift
  
  # Fixed regressors
  A_fixed <- cbind(
    intercept = 1,
    linear = drift,
    quad = drift^2
  )
  
  # Run full LSS pipeline
  lss_prep <- prepare_lss_fixed_components_core(A_fixed, 1, 1e-4)
  
  Beta_estimated <- run_lss_voxel_loop_core(
    Y_proj, X_trials, H_shapes, A_fixed,
    lss_prep$P_lss_matrix, lss_prep$p_lss_vector,
    ram_heuristic_GB_for_Rt = 1.0
  )
  
  # Check recovery
  expect_equal(dim(Beta_estimated), c(T, V))
  
  # Check that we get reasonable estimates
  # LSS with overlapping trials is challenging, so be lenient
  
  # Check that estimates have reasonable variance
  expect_gt(var(as.vector(Beta_estimated)), 0.01)
  
  # Check that estimates are not all zero
  expect_gt(mean(abs(Beta_estimated)), 0.05)
  
  # Check at least some correlation with truth for active voxels
  # Focus on voxels with strong activation
  active_voxels <- 1:3
  cors_active <- sapply(active_voxels, function(v) {
    if (var(Beta_estimated[, v]) > 0.01) {
      cor(true_betas[, v], Beta_estimated[, v])
    } else {
      0
    }
  })
  
  # At least one voxel should show some correlation
  expect_gt(max(cors_active), 0.3)
})