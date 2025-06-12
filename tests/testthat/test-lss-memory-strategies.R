# Test memory strategy implementations in consolidated core_lss.R
library(testthat)
library(manifoldhrf)

test_that("LSS memory strategies produce identical results", {
  # Create small test data
  set.seed(123)
  n <- 100  # timepoints
  V <- 10   # voxels  
  T_trials <- 5  # trials
  p <- 15   # HRF length
  m <- 3    # manifold dimensions
  
  # Create test data
  Y_proj <- matrix(rnorm(n * V), n, V)
  
  # Create trial design matrices
  X_trial_list <- lapply(1:T_trials, function(t) {
    X <- matrix(0, n, p)
    onset <- 5 + (t-1) * 18
    if (onset + p <= n) {
      X[onset:(onset + p - 1), ] <- diag(p)
    }
    X
  })
  
  # Create manifold components
  B_reconstructor <- matrix(rnorm(p * m), p, m)
  B_reconstructor <- qr.Q(qr(B_reconstructor))  # Orthonormalize
  
  Xi_smoothed <- matrix(rnorm(m * V), m, V)
  
  # Test all three strategies
  result_full <- run_lss_voxel_loop_core(
    Y_proj, X_trial_list, B_reconstructor, Xi_smoothed,
    memory_strategy = "full",
    verbose = FALSE
  )
  
  result_chunked <- run_lss_voxel_loop_core(
    Y_proj, X_trial_list, B_reconstructor, Xi_smoothed,
    memory_strategy = "chunked",
    chunk_size = 2,
    verbose = FALSE
  )
  
  result_streaming <- run_lss_voxel_loop_core(
    Y_proj, X_trial_list, B_reconstructor, Xi_smoothed,
    memory_strategy = "streaming",
    verbose = FALSE
  )
  
  # All strategies should produce nearly identical results
  # (small differences due to computation order are acceptable)
  expect_equal(result_full, result_chunked, tolerance = 1e-2)
  expect_equal(result_full, result_streaming, tolerance = 1e-2)
  
  # Check dimensions
  expect_equal(dim(result_full), c(T_trials, V))
})

test_that("Memory strategy auto-selection works correctly", {
  set.seed(456)
  
  # Small data - should select "full"
  small_result <- run_lss_voxel_loop_core(
    Y_proj_matrix = matrix(rnorm(50 * 10), 50, 10),
    X_trial_onset_list_of_matrices = lapply(1:3, function(i) matrix(rnorm(50 * 10), 50, 10)),
    B_reconstructor_matrix = matrix(rnorm(10 * 3), 10, 3),
    Xi_smoothed_allvox_matrix = matrix(rnorm(3 * 10), 3, 10),
    memory_strategy = "auto",
    ram_limit_GB = 4,
    verbose = TRUE
  )
  expect_equal(dim(small_result), c(3, 10))
  
  # Large data - should select "chunked" or "streaming"
  large_result <- run_lss_voxel_loop_core(
    Y_proj_matrix = matrix(rnorm(100 * 1000), 100, 1000),
    X_trial_onset_list_of_matrices = lapply(1:50, function(i) matrix(0, 100, 10)),
    B_reconstructor_matrix = matrix(rnorm(10 * 3), 10, 3),
    Xi_smoothed_allvox_matrix = matrix(rnorm(3 * 1000), 3, 1000),
    memory_strategy = "auto",
    ram_limit_GB = 0.1,  # Very low limit to force streaming
    verbose = TRUE
  )
  expect_equal(dim(large_result), c(50, 1000))
})

test_that("prepare_lss_fixed_components_core handles inputs correctly", {
  # No fixed regressors
  result_null <- prepare_lss_fixed_components_core(NULL)
  expect_null(result_null$P_lss)
  expect_false(result_null$has_intercept)
  
  # With fixed regressors including intercept
  n <- 100
  A_fixed <- cbind(
    1,  # intercept
    (1:n) / n,  # linear trend
    rnorm(n)  # noise regressor
  )
  
  result_fixed <- prepare_lss_fixed_components_core(A_fixed)
  expect_true(result_fixed$has_intercept)
  expect_null(result_fixed$P_lss)  # fmrilss handles internally
})

test_that("reconstruct_hrf_shapes_core performs correct matrix multiplication", {
  p <- 20
  m <- 4
  V <- 50
  
  B <- matrix(rnorm(p * m), p, m)
  Xi <- matrix(rnorm(m * V), m, V)
  
  H <- reconstruct_hrf_shapes_core(B, Xi)
  
  expect_equal(dim(H), c(p, V))
  expect_equal(H, B %*% Xi)
})

test_that("run_lss_for_voxel_core handles single voxel correctly", {
  set.seed(789)
  n <- 100
  T_trials <- 4
  p <- 15
  
  # Create single voxel data
  y_voxel <- rnorm(n)
  
  # Create trial matrices
  X_trials <- lapply(1:T_trials, function(t) {
    X <- matrix(0, n, p)
    onset <- 10 + (t-1) * 20
    if (onset + p <= n) {
      X[onset:(onset + p - 1), ] <- diag(p)
    }
    X
  })
  
  # Create HRF
  h_voxel <- dgamma(0:(p-1), shape = 5, rate = 1.5)
  h_voxel <- h_voxel / sum(h_voxel)
  
  # Run LSS
  betas <- run_lss_for_voxel_core(y_voxel, X_trials, h_voxel)
  
  expect_length(betas, T_trials)
  expect_type(betas, "double")
})

test_that("Validation functions catch errors correctly", {
  # Test validate_design_matrix_list
  expect_silent(validate_design_matrix_list(
    list(matrix(1:20, 10, 2), matrix(21:40, 10, 2)), 
    10
  ))
  
  expect_error(
    validate_design_matrix_list("not a list", 10),
    "Design matrices must be provided as a non-empty list"
  )
  
  # Test validate_hrf_shape_matrix
  expect_silent(validate_hrf_shape_matrix(
    matrix(1:50, 10, 5), 10, 5
  ))
  
  expect_error(
    validate_hrf_shape_matrix(matrix(1:50, 10, 5), 10, 6),
    "HRF shape matrix must have 6 columns, not 5"
  )
})

test_that("User-facing wrappers work correctly", {
  set.seed(111)
  n <- 80
  V <- 5
  T_trials <- 3
  p <- 12
  
  # Create test data
  Y <- matrix(rnorm(n * V), n, V)
  
  # Create trial list
  X_trials <- lapply(1:T_trials, function(t) {
    X <- matrix(0, n, p)
    onset <- 5 + (t-1) * 25
    if (onset + p <= n) {
      X[onset:(onset + p - 1), ] <- diag(p)
    }
    X
  })
  
  # Create HRF matrix
  H <- matrix(0, p, V)
  for (v in 1:V) {
    h <- dgamma(0:(p-1), shape = 4 + v/2, rate = 1.5)
    H[, v] <- h / sum(h)
  }
  
  # Test run_lss_voxel_loop wrapper
  result <- run_lss_voxel_loop(Y, X_trials, H, 
                               memory_strategy = "full",
                               verbose = FALSE)
  
  expect_equal(dim(result), c(T_trials, V))
  
  # Test single voxel wrapper
  single_result <- run_lss_for_voxel(Y[,1], X_trials, H[,1])
  expect_length(single_result, T_trials)
})