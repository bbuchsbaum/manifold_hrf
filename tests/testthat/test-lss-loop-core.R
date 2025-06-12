# Tests for run_lss_voxel_loop_core

test_that("run_lss_voxel_loop_core matches single voxel implementation", {
  set.seed(123)
  n <- 40
  p <- 8
  V <- 3
  T_trials <- 5
  m <- 3  # manifold dimensions

  # Create manifold components
  B_reconstructor <- matrix(rnorm(p * m), p, m)
  B_reconstructor <- qr.Q(qr(B_reconstructor))  # Orthonormalize
  Xi_smoothed <- matrix(rnorm(m * V), m, V)
  
  # HRF shapes for each voxel (reconstructed from manifold)
  H_shapes <- B_reconstructor %*% Xi_smoothed

  # Trial design matrices
  X_trials <- lapply(seq_len(T_trials), function(t) {
    X <- matrix(0, n, p)
    onset <- sample(1:(n - p), 1)
    for (j in seq_len(p)) {
      X[onset + j - 1, j] <- 1
    }
    X
  })

  # True trial amplitudes
  Beta_true <- matrix(rnorm(T_trials * V), T_trials, V)

  # Confounds and projection
  A_fixed <- cbind(1, rnorm(n))
  # Manual projection matrix calculation for data generation
  AtA <- crossprod(A_fixed)
  AtA_reg <- AtA + 1e-6 * diag(ncol(A_fixed))
  P_conf <- diag(n) - A_fixed %*% solve(AtA_reg) %*% t(A_fixed)

  # Generate projected data
  Y_clean <- matrix(0, n, V)
  for (v in seq_len(V)) {
    for (t in seq_len(T_trials)) {
      Y_clean[, v] <- Y_clean[, v] + X_trials[[t]] %*% (H_shapes[, v] * Beta_true[t, v])
    }
  }
  Y_proj <- P_conf %*% (Y_clean + matrix(rnorm(n * V, sd = 0.05), n, V))

  # Use the current API for run_lss_voxel_loop_core
  Beta_core <- run_lss_voxel_loop_core(
    Y_proj_matrix = Y_proj,
    X_trial_onset_list_of_matrices = X_trials,
    B_reconstructor_matrix = B_reconstructor,
    Xi_smoothed_allvox_matrix = Xi_smoothed,
    A_lss_fixed_matrix = A_fixed,
    memory_strategy = "full",
    verbose = FALSE
  )

  # Manual implementation using the same components
  Beta_manual <- matrix(0, T_trials, V)
  for (v in seq_len(V)) {
    # Use run_lss_for_voxel_core for consistency (same underlying implementation)
    Beta_manual[, v] <- run_lss_for_voxel_core(
      Y_proj_voxel_vector = Y_proj[, v],
      X_trial_onset_list_of_matrices = X_trials,
      h_voxel_shape_vector = H_shapes[, v],
      A_lss_fixed_matrix = A_fixed
    )
  }

  expect_equal(Beta_core, Beta_manual, tolerance = 1e-5)
})

test_that("run_lss_voxel_loop_core works with streaming strategy", {
  set.seed(123)
  n <- 40
  p <- 8
  V <- 3
  T_trials <- 5
  m <- 3  # manifold dimensions

  # Create manifold components
  B_reconstructor <- matrix(rnorm(p * m), p, m)
  B_reconstructor <- qr.Q(qr(B_reconstructor))  # Orthonormalize
  Xi_smoothed <- matrix(rnorm(m * V), m, V)
  
  # HRF shapes for each voxel (reconstructed from manifold)
  H_shapes <- B_reconstructor %*% Xi_smoothed
  
  X_trials <- lapply(seq_len(T_trials), function(t) {
    X <- matrix(0, n, p)
    onset <- sample(1:(n - p), 1)
    for (j in seq_len(p)) {
      X[onset + j - 1, j] <- 1
    }
    X
  })

  Beta_true <- matrix(rnorm(T_trials * V), T_trials, V)

  A_fixed <- cbind(1, rnorm(n))
  # Manual projection matrix calculation
  AtA <- crossprod(A_fixed)
  AtA_reg <- AtA + 1e-6 * diag(ncol(A_fixed))
  P_conf <- diag(n) - A_fixed %*% solve(AtA_reg) %*% t(A_fixed)

  Y_clean <- matrix(0, n, V)
  for (v in seq_len(V)) {
    for (t in seq_len(T_trials)) {
      Y_clean[, v] <- Y_clean[, v] + X_trials[[t]] %*% (H_shapes[, v] * Beta_true[t, v])
    }
  }
  Y_proj <- P_conf %*% (Y_clean + matrix(rnorm(n * V, sd = 0.05), n, V))

  # Test streaming strategy
  Beta_streaming <- run_lss_voxel_loop_core(
    Y_proj_matrix = Y_proj,
    X_trial_onset_list_of_matrices = X_trials,
    B_reconstructor_matrix = B_reconstructor,
    Xi_smoothed_allvox_matrix = Xi_smoothed,
    A_lss_fixed_matrix = A_fixed,
    memory_strategy = "streaming",
    verbose = FALSE
  )

  # Test full strategy for comparison
  Beta_full <- run_lss_voxel_loop_core(
    Y_proj_matrix = Y_proj,
    X_trial_onset_list_of_matrices = X_trials,
    B_reconstructor_matrix = B_reconstructor,
    Xi_smoothed_allvox_matrix = Xi_smoothed,
    A_lss_fixed_matrix = A_fixed,
    memory_strategy = "full",
    verbose = FALSE
  )

  # Different strategies should give the same result
  expect_equal(Beta_streaming, Beta_full, tolerance = 1e-5)
})
