# Tests for Core LSS Functions with fmrilss backend

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
    A_fixed_regressors_matrix = A_fixed, 
    lambda_ridge_A = 1e-6
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("P_lss", "has_intercept"))
  
  # Check that it detected the intercept
  expect_true(result$has_intercept)
  
  # Check that P_lss is NULL (fmrilss handles projection internally)
  expect_null(result$P_lss)
})

test_that("reconstruct_hrf_shapes_core works correctly", {
  set.seed(789)
  
  # HRF library dimension
  p <- 20  # HRF samples
  m <- 3   # manifold dimensions
  V <- 5   # voxels
  
  # Create test reconstructor and coordinates
  B_reconstructor <- matrix(rnorm(p * m), p, m)
  Xi_smoothed <- matrix(rnorm(m * V), m, V)
  
  # Reconstruct HRF shapes
  H_shapes <- reconstruct_hrf_shapes_core(B_reconstructor, Xi_smoothed)
  
  # Check dimensions
  expect_equal(dim(H_shapes), c(p, V))
  
  # Check that reconstruction is simply matrix multiplication
  expected <- B_reconstructor %*% Xi_smoothed
  expect_equal(H_shapes, expected)
})

test_that("run_lss_for_voxel works correctly", {
  set.seed(101)
  
  # Test parameters
  n <- 60
  p <- 10
  T_trials <- 4
  
  # Create trial matrices
  X_trials <- lapply(seq_len(T_trials), function(t) {
    X <- matrix(0, n, p)
    onset <- 5 + (t-1) * (p + 2)
    if (onset + p - 1 <= n) {
      X[onset:(onset + p - 1), ] <- diag(p)
    }
    X
  })
  
  # HRF shape
  h <- exp(-(0:(p-1))/3)
  h <- h / sum(h)
  
  # Generate data
  true_betas <- rnorm(T_trials)
  y <- rnorm(n, sd = 0.1)
  for (t in seq_len(T_trials)) {
    y <- y + X_trials[[t]] %*% h * true_betas[t]
  }
  
  # Run LSS
  result <- run_lss_for_voxel(
    y_voxel = y,
    X_trial_list = X_trials,
    h_voxel = h,
    TR = 2
  )
  
  # Check output (now returns a vector directly)
  expect_type(result, "double")
  expect_length(result, T_trials)
  expect_true(all(is.finite(result)))
  
  # Estimates should be reasonably close to true values
  # (won't be exact due to noise and intercept-only model)
  expect_lt(mean(abs(result - true_betas)), 0.5)
})

test_that("run_lss_voxel_loop_core handles multiple voxels", {
  set.seed(202)
  
  n <- 40
  p <- 8
  V <- 3
  T_trials <- 2
  
  # Create data
  Y_proj <- matrix(rnorm(n * V), n, V)
  H_shapes <- matrix(rnorm(p * V), p, V)
  
  X_trials <- lapply(seq_len(T_trials), function(t) {
    X <- matrix(0, n, p)
    onset <- 5 + (t-1) * (p + 5)
    if (onset + p - 1 <= n) {
      X[onset:(onset + p - 1), ] <- diag(p)
    }
    X
  })
  
  # Fixed regressors (not used in current implementation)
  A_fixed <- cbind(1, rnorm(n))
  lss_prep <- prepare_lss_fixed_components_core(
    A_fixed_regressors_matrix = A_fixed, 
    lambda_ridge_A = 1e-6
  )
  
  # Create mock B_reconstructor and Xi_smoothed for new API
  m <- 3  # manifold dimensions
  B_reconstructor <- matrix(rnorm(p * m), p, m)
  Xi_smoothed <- matrix(rnorm(m * V), m, V)
  
  # Run voxel loop
  Beta <- run_lss_voxel_loop_core(
    Y_proj_matrix = Y_proj,
    X_trial_onset_list_of_matrices = X_trials,
    B_reconstructor_matrix = B_reconstructor,
    Xi_smoothed_allvox_matrix = Xi_smoothed,
    A_lss_fixed_matrix = A_fixed,
    memory_strategy = "auto",
    n_cores = 1,
    progress = FALSE,
    verbose = FALSE
  )
  
  # Check output
  expect_equal(dim(Beta), c(T_trials, V))
  expect_true(all(is.finite(Beta)))
})

test_that("fmrilss backend handles edge cases", {
  set.seed(303)
  
  # Edge case: single trial
  n <- 30
  p <- 5
  
  X_single <- list(matrix(rnorm(n * p), n, p))
  h <- rnorm(p)
  y <- rnorm(n)
  
  result <- run_lss_for_voxel(
    y_voxel = y,
    X_trial_list = X_single,
    h_voxel = h
  )
  
  expect_length(result, 1)
  expect_true(is.finite(result[1]))
  
  # Edge case: trials at boundaries
  X_boundary <- list(
    matrix(c(rep(1, p), rep(0, n-p)), n, 1),  # Trial at start
    matrix(c(rep(0, n-p), rep(1, p)), n, 1)   # Trial at end
  )
  
  result_boundary <- run_lss_for_voxel(
    y_voxel = y,
    X_trial_list = X_boundary,
    h_voxel = rep(1, 1)
  )
  
  expect_length(result_boundary, 2)
})

test_that("validate_design_matrix_list catches errors", {
  n <- 50
  
  # Valid list
  X_valid <- list(
    matrix(1:50, 50, 1),
    matrix(51:100, 50, 1)
  )
  expect_silent(validate_design_matrix_list(X_valid, n))
  
  # Not a list
  expect_error(
    validate_design_matrix_list(matrix(1:50, 50, 1), n),
    "Design matrices must be provided as a non-empty list"
  )
  
  # Empty list
  expect_error(
    validate_design_matrix_list(list(), n),
    "Design matrices must be provided as a non-empty list"
  )
  
  # Wrong dimensions
  X_wrong <- list(
    matrix(1:40, 40, 1),  # Wrong number of rows
    matrix(1:50, 50, 1)
  )
  expect_error(
    validate_design_matrix_list(X_wrong, n),
    "Design matrix 1 has 40 rows but expected 50"
  )
})