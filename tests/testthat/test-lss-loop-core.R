# Tests for run_lss_voxel_loop_core

test_that("run_lss_voxel_loop_core matches single voxel implementation", {
  set.seed(123)
  n <- 40
  p <- 8
  V <- 3
  T_trials <- 5

  # HRF shapes for each voxel
  H_shapes <- matrix(rnorm(p * V), p, V)

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
  lss_prep <- prepare_lss_fixed_components_core(A_fixed, 1, 1e-6)
  P_conf <- prepare_projection_matrix(A_fixed, 1e-6)

  # Generate projected data
  Y_clean <- matrix(0, n, V)
  for (v in seq_len(V)) {
    for (t in seq_len(T_trials)) {
      Y_clean[, v] <- Y_clean[, v] + X_trials[[t]] %*% (H_shapes[, v] * Beta_true[t, v])
    }
  }
  Y_proj <- P_conf %*% (Y_clean + matrix(rnorm(n * V, sd = 0.05), n, V))

  Beta_core <- manifoldhrf:::run_lss_voxel_loop_core(
    Y_proj_matrix = Y_proj,
    X_trial_onset_list_of_matrices = X_trials,
    H_shapes_allvox_matrix = H_shapes,
    A_lss_fixed_matrix = A_fixed,
    P_lss_matrix = lss_prep$P_lss_matrix,
    p_lss_vector = lss_prep$p_lss_vector,
    n_jobs = 1,
    ram_heuristic_GB_for_Rt = 1.0
  )

  Beta_manual <- matrix(0, T_trials, V)
  for (v in seq_len(V)) {
    Beta_manual[, v] <- manifoldhrf::run_lss_for_voxel_corrected_full(
      Y_proj_voxel_vector = Y_proj[, v],
      X_trial_onset_list_of_matrices = X_trials,
      H_shape_voxel_vector = H_shapes[, v],
      A_lss_fixed_matrix = A_fixed,
      P_lss_matrix = lss_prep$P_lss_matrix,
      p_lss_vector = lss_prep$p_lss_vector
    )
  }

  expect_equal(Beta_core, Beta_manual, tolerance = 1e-8)
})

test_that("run_lss_voxel_loop_core works without precomputation", {
  set.seed(123)
  n <- 40
  p <- 8
  V <- 3
  T_trials <- 5

  H_shapes <- matrix(rnorm(p * V), p, V)
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
  lss_prep <- prepare_lss_fixed_components_core(A_fixed, 1, 1e-6)
  P_conf <- prepare_projection_matrix(A_fixed, 1e-6)

  Y_clean <- matrix(0, n, V)
  for (v in seq_len(V)) {
    for (t in seq_len(T_trials)) {
      Y_clean[, v] <- Y_clean[, v] + X_trials[[t]] %*% (H_shapes[, v] * Beta_true[t, v])
    }
  }
  Y_proj <- P_conf %*% (Y_clean + matrix(rnorm(n * V, sd = 0.05), n, V))

  Beta_core <- manifoldhrf:::run_lss_voxel_loop_core(
    Y_proj_matrix = Y_proj,
    X_trial_onset_list_of_matrices = X_trials,
    H_shapes_allvox_matrix = H_shapes,
    A_lss_fixed_matrix = A_fixed,
    P_lss_matrix = lss_prep$P_lss_matrix,
    p_lss_vector = lss_prep$p_lss_vector,
    n_jobs = 1,
    ram_heuristic_GB_for_Rt = 0
  )

  Beta_manual <- matrix(0, T_trials, V)
  for (v in seq_len(V)) {
    Beta_manual[, v] <- run_lss_woodbury_corrected(
      Y_proj_voxel_vector = Y_proj[, v],
      X_trial_onset_list_of_matrices = X_trials,
      H_shape_voxel_vector = H_shapes[, v],
      P_confound = P_conf,
      lambda_ridge = 1e-6
    )
  }

  expect_equal(Beta_core, Beta_manual, tolerance = 1e-8)
})
