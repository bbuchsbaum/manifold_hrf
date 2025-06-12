# Test LS-A vs LS-S equivalence for T=2

test_that("LS-A and LS-S are equivalent for T=2", {
  set.seed(42)
  
  # Create synthetic data
  n <- 100  # timepoints
  T_trials <- 2  # exactly 2 trials for equivalence
  p <- 3    # HRF basis functions
  
  # Trial design matrices
  X_trial_list <- list()
  for (t in 1:T_trials) {
    X_trial_list[[t]] <- matrix(rnorm(n * p), n, p)
  }
  
  # HRF shape
  h_voxel <- rnorm(p)
  
  # Response data
  y_voxel <- rnorm(n)
  
  # No confounds case - use direct LS-A computation
  # Build C matrix for LS-A
  C_lsa <- matrix(0, n, T_trials)
  for (t in seq_len(T_trials)) {
    C_lsa[, t] <- X_trial_list[[t]] %*% h_voxel
  }
  
  # LS-A: solve normal equations directly (since no confounds)
  lsa_result <- solve(crossprod(C_lsa), crossprod(C_lsa, y_voxel))
  
  lss_result <- run_lss_for_voxel(
    y_voxel = y_voxel,
    X_trial_list = X_trial_list,
    h_voxel = h_voxel
  )
  
  expect_equal(
    as.vector(lsa_result),
    as.vector(lss_result),
    tolerance = 1e-2,
    info = "LS-A and LS-S should be equivalent for T=2"
  )
})

test_that("LS-A and LS-S are equivalent for T=2 with confounds", {
  set.seed(123)
  
  # Create synthetic data
  n <- 80
  T_trials <- 2
  p <- 2
  q <- 3  # confounds
  
  # Trial design matrices
  X_trial_list <- list()
  for (t in 1:T_trials) {
    X_trial_list[[t]] <- matrix(rnorm(n * p), n, p)
  }
  
  # HRF shape and confounds
  h_voxel <- rnorm(p)
  Z_confounds <- matrix(rnorm(n * q), n, q)
  y_voxel <- rnorm(n)
  
  # Build C matrix
  C <- matrix(0, n, T_trials)
  for (t in seq_len(T_trials)) {
    C[, t] <- X_trial_list[[t]] %*% h_voxel
  }
  
  # LS-S result using fmrilss directly with confounds
  lss_result <- fmrilss::lss(
    Y = matrix(y_voxel, ncol = 1),
    X = C,
    Z = Z_confounds,
    method = "r_optimized"
  )
  
  # LS-A result: For T=2, LSS and LSA should be equivalent
  # Combine confounds with trial regressors for full model
  X_full <- cbind(C, Z_confounds)
  lsa_betas_full <- solve(crossprod(X_full), crossprod(X_full, y_voxel))
  lsa_betas <- lsa_betas_full[1:T_trials]
  
  expect_equal(
    as.vector(lsa_betas),
    as.vector(lss_result),
    tolerance = 1e-2
  )
})

test_that("Memory usage is reasonable for LS-S", {
  skip_if_not_installed("pryr")
  
  n <- 300  # typical fMRI session
  T_trials <- 50
  p <- 3
  
  # Create data
  X_trial_list <- list()
  for (t in 1:T_trials) {
    X_trial_list[[t]] <- matrix(rnorm(n * p), n, p)
  }
  h_voxel <- rnorm(p)
  y_voxel <- rnorm(n)
  
  # Test memory usage
  mem_change <- pryr::mem_change({
    result <- run_lss_for_voxel(
      y_voxel = y_voxel,
      X_trial_list = X_trial_list,
      h_voxel = h_voxel
    )
  })
  
  # Should not allocate large n×n matrices
  expect_lt(
    abs(mem_change), 
    1e6  # 1 MB
  )
}) 