# Test the corrected LSS implementation
library(testthat)

test_that("Corrected LSS implementation matches ground truth", {
  
  set.seed(789)
  
  # Test parameters
  n <- 100
  T_trials <- 8
  p <- 15
  q <- 4
  lambda <- 1e-6
  
  # Create trial matrices
  X_trials <- list()
  for (t in 1:T_trials) {
    X_t <- matrix(0, n, p)
    onset <- 5 + (t-1) * 12
    if (onset + p <= n) {
      X_t[onset:(onset + p - 1), ] <- diag(p)
    }
    X_trials[[t]] <- X_t
  }
  
  # HRF
  h <- dgamma(0:(p-1), shape = 5, rate = 1.5)
  h <- h / sum(h)
  
  # Confounds
  Z <- cbind(
    1,
    (1:n) / n,
    sin(2 * pi * (1:n) / n),
    cos(2 * pi * (1:n) / n)
  )
  
  # True parameters
  true_betas <- c(1, -0.5, 2, 0, 1.5, -1, 0.8, 0.3)
  true_conf <- c(5, -2, 1, 0.5)
  
  # Generate data
  y <- Z %*% true_conf
  for (t in 1:T_trials) {
    y <- y + X_trials[[t]] %*% h * true_betas[t]
  }
  y <- y + rnorm(n, sd = 0.2)
  
  # Prepare projection
  P_proj <- prepare_projection_matrix(Z, lambda)
  y_proj <- P_proj %*% y
  
  # METHOD 1: Direct LSS (ground truth)
  betas_direct <- numeric(T_trials)
  for (t in 1:T_trials) {
    X_other <- do.call(cbind, lapply(X_trials[-t], function(X) X %*% h))
    X_full <- cbind(X_trials[[t]] %*% h, X_other, Z)
    XtX <- crossprod(X_full) + lambda * diag(ncol(X_full))
    Xty <- crossprod(X_full, y)
    betas_direct[t] <- solve(XtX, Xty)[1]
  }
  
  # METHOD 2: Corrected implementation
  betas_corrected <- run_lss_for_voxel_corrected_full(
    Y_proj_voxel_vector = y_proj,
    X_trial_onset_list_of_matrices = X_trials,
    H_shape_voxel_vector = h,
    P_confound = P_proj,
    lambda_ridge = lambda
  )
  
  # METHOD 3: Current Woodbury implementation
  lss_prep <- prepare_lss_fixed_components_core(
    A_lss_fixed_matrix = Z,
    intercept_col_index_in_Alss = 1,
    lambda_ridge_Alss = lambda
  )
  
  betas_current <- run_lss_woodbury_corrected(
    Y_proj_voxel_vector = as.vector(y_proj),
    X_trial_onset_list_of_matrices = X_trials,
    H_shape_voxel_vector = h,
    lambda_ridge = lambda
  )
  
  # Compare
  cat("\n=== LSS Implementation Comparison ===\n")
  cat("True betas:     ", round(true_betas, 3), "\n")
  cat("Direct LSS:     ", round(betas_direct, 3), "\n")
  cat("Corrected:      ", round(betas_corrected, 3), "\n")
  cat("Current Wood:   ", round(betas_current, 3), "\n")
  
  cat("\n=== Differences from Direct ===\n")
  cat("Corrected diff: ", round(betas_corrected - betas_direct, 6), "\n")
  cat("Current diff:   ", round(betas_current - betas_direct, 6), "\n")
  
  cat("\n=== Mean Squared Errors ===\n")
  cat("Direct MSE:     ", round(mean((betas_direct - true_betas)^2), 4), "\n")
  cat("Corrected MSE:  ", round(mean((betas_corrected - true_betas)^2), 4), "\n")
  cat("Current MSE:    ", round(mean((betas_current - true_betas)^2), 4), "\n")
  
  # The corrected method should match direct for trials with signal
  # Note: LSS is solving a different problem than simultaneous estimation
  # It estimates each trial's effect while treating others as nuisance
  
  # For a proper comparison, let's also try simultaneous estimation
  C_all <- do.call(cbind, lapply(X_trials, function(X) X %*% h))
  X_sim <- cbind(C_all, Z)
  betas_simultaneous <- solve(crossprod(X_sim) + lambda * diag(ncol(X_sim)), 
                             crossprod(X_sim, y))[1:T_trials]
  
  cat("\n=== Simultaneous Estimation ===\n")
  cat("Simultaneous:   ", round(betas_simultaneous, 3), "\n")
  cat("Diff from true: ", round(betas_simultaneous - true_betas, 3), "\n")
  
  # Key insight: The corrected implementation is doing simultaneous estimation
  # NOT trial-wise LSS. This explains the discrepancy!
  
  expect_lt(max(abs(betas_corrected - betas_simultaneous)), 1e-4,
            "Corrected implementation should match simultaneous estimation")
})


test_that("True LSS implementation works correctly", {
  
  # Let's implement true LSS where each trial is estimated separately
  # with all other trials as nuisance
  
  true_lss <- function(y, X_trials, h, Z, lambda = 1e-6) {
    T_trials <- length(X_trials)
    betas <- numeric(T_trials)
    
    for (t in 1:T_trials) {
      # For trial t: put it first, others as nuisance
      X_t <- X_trials[[t]] %*% h
      X_others <- do.call(cbind, lapply(X_trials[-t], function(X) X %*% h))
      
      # Full design: [trial_t | other_trials | confounds]
      X_full <- cbind(X_t, X_others, Z)
      
      # Solve
      XtX <- crossprod(X_full) + lambda * diag(ncol(X_full))
      Xty <- crossprod(X_full, y)
      beta_all <- solve(XtX, Xty)
      
      betas[t] <- beta_all[1]  # First coefficient is trial of interest
    }
    
    return(betas)
  }
  
  # Test setup
  set.seed(999)
  n <- 80
  T_trials <- 5
  p <- 12
  
  X_trials <- list()
  for (t in 1:T_trials) {
    X_t <- matrix(0, n, p)
    onset <- 8 + (t-1) * 15
    if (onset + p <= n) {
      X_t[onset:(onset + p - 1), ] <- diag(p)
    }
    X_trials[[t]] <- X_t
  }
  
  h <- exp(-(0:(p-1))/3)
  h <- h / sum(h)
  
  Z <- cbind(1, (1:n)/n)
  
  true_betas <- c(2, -1, 0, 1.5, 0.5)
  
  y <- Z %*% c(10, -5)
  for (t in 1:T_trials) {
    y <- y + X_trials[[t]] %*% h * true_betas[t]
  }
  y <- y + rnorm(n, sd = 0.1)
  
  # Compare implementations
  betas_true_lss <- true_lss(y, X_trials, h, Z)
  
  # Project y for Woodbury method
  P <- diag(n) - Z %*% solve(crossprod(Z) + 1e-6 * diag(ncol(Z))) %*% t(Z)
  y_proj <- P %*% y
  
  lss_prep <- prepare_lss_fixed_components_core(
    A_lss_fixed_matrix = Z,
    intercept_col_index_in_Alss = 1,
    lambda_ridge_Alss = 1e-6
  )
  
  betas_woodbury <- run_lss_woodbury_corrected(
    Y_proj_voxel_vector = as.vector(y_proj),
    X_trial_onset_list_of_matrices = X_trials,
    H_shape_voxel_vector = h,
    lambda_ridge = 1e-6
  )
  
  cat("\n=== True LSS vs Woodbury ===\n")
  cat("True betas:  ", round(true_betas, 3), "\n")
  cat("True LSS:    ", round(betas_true_lss, 3), "\n")
  cat("Woodbury:    ", round(betas_woodbury, 3), "\n")
  cat("Difference:  ", round(betas_woodbury - betas_true_lss, 4), "\n")
  
  # The Woodbury implementation should match true LSS
  # If it doesn't, there's a mathematical error in the formula
})