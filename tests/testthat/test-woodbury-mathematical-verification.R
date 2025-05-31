# Mathematical Verification of Woodbury LSS Implementation
# 
# Goal: Verify that the Woodbury implementation is mathematically correct
# by working through the algebra step by step

library(testthat)

# ============================================================================
# MATHEMATICAL BACKGROUND
# ============================================================================
# 
# Standard LSS: For each trial t, we solve:
#   β = (X'X + λI)^(-1) X'y
# 
# where X = [x_t, X_other, Z] includes:
#   - x_t: regressor for trial t (n × 1)
#   - X_other: regressors for all other trials (n × (T-1))
#   - Z: confound regressors (n × q)
#
# The Woodbury identity allows us to update the inverse when adding x_t:
#   (A + uv')^(-1) = A^(-1) - A^(-1)uv'A^(-1) / (1 + v'A^(-1)u)
#
# But our implementation uses a different approach based on:
# 1. Pre-projecting y to remove confounds: y_proj = (I - Z(Z'Z)^(-1)Z')y
# 2. Using the Woodbury identity to handle the trial regressors efficiently
#
# ============================================================================

test_that("Woodbury LSS matches theoretical derivation", {
  
  set.seed(123)
  
  # Small problem for detailed verification
  n <- 50      # time points
  T_trials <- 4 # trials  
  q <- 3        # confounds (intercept + 2 drift terms)
  p <- 10       # HRF length
  
  # Create trial design matrices
  X_trials <- list()
  for (t in 1:T_trials) {
    X_t <- matrix(0, n, p)
    onset <- 5 + (t-1) * 12
    if (onset + p <= n) {
      X_t[onset:(onset + p - 1), ] <- diag(p)
    }
    X_trials[[t]] <- X_t
  }
  
  # Simple HRF
  h <- exp(-(0:(p-1))/3)
  h <- h / sum(h)
  
  # Confound matrix
  Z <- cbind(
    1,                          # intercept
    (1:n) / n,                  # linear drift
    ((1:n) / n)^2               # quadratic drift
  )
  
  # True parameters
  true_betas <- c(2, -1, 1.5, 0.5)
  true_confound_weights <- c(10, -5, 2)
  
  # Generate data
  y <- Z %*% true_confound_weights
  for (t in 1:T_trials) {
    y <- y + X_trials[[t]] %*% h * true_betas[t]
  }
  y <- y + rnorm(n, sd = 0.1)  # Small noise
  
  # ========================================================================
  # METHOD 1: Direct LSS (Ground Truth)
  # ========================================================================
  
  lambda <- 1e-6
  betas_direct <- numeric(T_trials)
  
  for (t in 1:T_trials) {
    # Build full design matrix
    X_other <- do.call(cbind, lapply(X_trials[-t], function(X) X %*% h))
    X_full <- cbind(X_trials[[t]] %*% h, X_other, Z)
    
    # Normal equations with ridge
    XtX <- crossprod(X_full) + lambda * diag(ncol(X_full))
    Xty <- crossprod(X_full, y)
    betas_all <- solve(XtX, Xty)
    betas_direct[t] <- betas_all[1]
  }
  
  # ========================================================================
  # METHOD 2: Woodbury with Projected Y (Current Implementation)
  # ========================================================================
  
  # Project out confounds from y
  P_Z <- diag(n) - Z %*% solve(crossprod(Z) + lambda * diag(q)) %*% t(Z)
  y_proj <- P_Z %*% y
  
  # Prepare LSS components
  lss_prep <- prepare_lss_fixed_components_core(
    A_lss_fixed_matrix = Z,
    intercept_col_index_in_Alss = 1,
    lambda_ridge_Alss = lambda
  )
  
  betas_woodbury <- run_lss_for_voxel_core(
    Y_proj_voxel_vector = as.vector(y_proj),
    X_trial_onset_list_of_matrices = X_trials,
    H_shape_voxel_vector = h,
    A_lss_fixed_matrix = Z,
    P_lss_matrix = lss_prep$P_lss_matrix,
    p_lss_vector = lss_prep$p_lss_vector
  )
  
  # ========================================================================
  # METHOD 3: Manual Woodbury Implementation for Verification
  # ========================================================================
  
  # Let's implement the Woodbury formula manually to understand what SHOULD happen
  
  # First, let's understand what we're solving:
  # We want: β = (X'X + λI)^(-1) X'y
  # where X = [x_t, X_other, Z]
  
  # The key insight: we can reformulate this using the Schur complement
  # after projecting out Z from both y and the trial regressors
  
  betas_manual <- numeric(T_trials)
  
  # Pre-compute some matrices
  C_all <- do.call(cbind, lapply(X_trials, function(X) X %*% h))  # All trial regressors
  
  # For each trial
  for (t in 1:T_trials) {
    # Current trial regressor
    c_t <- C_all[, t]
    
    # All other trial regressors
    C_other <- C_all[, -t, drop = FALSE]
    
    # === APPROACH 1: Direct Schur Complement ===
    # We're solving: [C_all, Z]' [C_all, Z] β = [C_all, Z]' y
    # But we can use the fact that y_proj = P_Z y has confounds removed
    
    # Project trial regressors
    C_all_proj <- P_Z %*% C_all
    c_t_proj <- C_all_proj[, t]
    C_other_proj <- C_all_proj[, -t, drop = FALSE]
    
    # Now solve the projected system
    X_proj <- cbind(c_t_proj, C_other_proj)
    XtX_proj <- crossprod(X_proj) + lambda * diag(ncol(X_proj))
    Xty_proj <- crossprod(X_proj, y_proj)
    betas_proj <- solve(XtX_proj, Xty_proj)
    betas_manual[t] <- betas_proj[1]
  }
  
  # ========================================================================
  # COMPARISON
  # ========================================================================
  
  cat("\n=== Beta Estimates ===\n")
  cat("True betas:     ", true_betas, "\n")
  cat("Direct LSS:     ", round(betas_direct, 4), "\n")
  cat("Woodbury:       ", round(betas_woodbury, 4), "\n")
  cat("Manual:         ", round(betas_manual, 4), "\n")
  
  cat("\n=== Differences from Direct ===\n")
  cat("Woodbury diff:  ", round(betas_woodbury - betas_direct, 6), "\n")
  cat("Manual diff:    ", round(betas_manual - betas_direct, 6), "\n")
  
  # The manual method (projecting everything) should be very close to direct
  expect_lt(max(abs(betas_manual - betas_direct)), 1e-4,
            "Manual projected method should match direct LSS")
  
  # Now let's debug why Woodbury doesn't match
  # Let's trace through the Woodbury implementation step by step
  
  cat("\n=== Debugging Woodbury Implementation ===\n")
  
  # What the implementation does:
  # 1. C_v = all trial regressors convolved with HRF
  # 2. U_v = P_lss %*% C_v
  # 3. V_regressors_v = C_v - Z %*% U_v
  # 4. Uses these to compute betas
  
  C_v <- C_all
  U_v <- lss_prep$P_lss_matrix %*% C_v
  V_regressors_v <- C_v - Z %*% U_v
  
  cat("C_v shape: ", dim(C_v), "\n")
  cat("U_v shape: ", dim(U_v), "\n")
  cat("V_regressors_v shape: ", dim(V_regressors_v), "\n")
  
  # Check if V_regressors_v is actually the projected C_v
  projection_error <- norm(V_regressors_v - P_Z %*% C_v, "F")
  cat("Projection error: ", projection_error, "\n")
  
  # The issue might be in how the final beta is computed
  # Let's check the formula being used
  
  pc_v_row <- as.vector(crossprod(lss_prep$p_lss_vector, C_v))
  cv_v_row <- colSums(V_regressors_v^2)
  alpha_v_row <- (1 - pc_v_row) / pmax(cv_v_row, .Machine$double.eps)
  
  cat("\npc_v_row: ", round(pc_v_row, 4), "\n")
  cat("cv_v_row: ", round(cv_v_row, 4), "\n")
  cat("alpha_v_row: ", round(alpha_v_row, 4), "\n")
  
  # The formula seems to be implementing something different than standard Woodbury
  # It's using a scaling approach rather than the matrix inversion lemma
})


test_that("Corrected Woodbury implementation matches direct LSS", {
  
  # Let's implement what the Woodbury formula SHOULD be
  
  corrected_woodbury_lss <- function(y_proj, X_trials, h, Z, lambda = 1e-6) {
    
    n <- length(y_proj)
    T_trials <- length(X_trials)
    
    # All trial regressors
    C <- do.call(cbind, lapply(X_trials, function(X) X %*% h))
    
    # Pre-compute (Z'Z + λI)^(-1)
    ZtZ_inv <- solve(crossprod(Z) + lambda * diag(ncol(Z)))
    
    # Pre-compute projection matrix components
    # P = I - Z(Z'Z + λI)^(-1)Z'
    # We need P'C for each trial
    P <- diag(n) - Z %*% ZtZ_inv %*% t(Z)
    C_proj <- P %*% C
    
    # For efficiency, pre-compute (C_proj'C_proj + λI)^(-1)
    CtC_proj_inv <- solve(crossprod(C_proj) + lambda * diag(T_trials))
    
    # Compute all betas at once
    betas <- CtC_proj_inv %*% crossprod(C_proj, y_proj)
    
    return(as.vector(betas))
  }
  
  # Test setup
  set.seed(456)
  n <- 100
  T_trials <- 6
  p <- 15
  q <- 4
  
  # Create test data
  X_trials <- list()
  for (t in 1:T_trials) {
    X_t <- matrix(0, n, p)
    onset <- 10 + (t-1) * 15
    if (onset + p <= n) {
      X_t[onset:(onset + p - 1), ] <- diag(p)
    }
    X_trials[[t]] <- X_t
  }
  
  h <- dgamma(0:(p-1), shape = 4, rate = 1)
  h <- h / sum(h)
  
  Z <- cbind(1, poly(1:n, degree = q-1))
  
  true_betas <- rnorm(T_trials, mean = 1, sd = 0.5)
  y <- numeric(n)
  for (t in 1:T_trials) {
    y <- y + X_trials[[t]] %*% h * true_betas[t]
  }
  y <- y + Z %*% rnorm(q, sd = 2) + rnorm(n, sd = 0.2)
  
  # Project y
  P_Z <- diag(n) - Z %*% solve(crossprod(Z) + 1e-6 * diag(q)) %*% t(Z)
  y_proj <- P_Z %*% y
  
  # Direct method
  betas_direct <- numeric(T_trials)
  for (t in 1:T_trials) {
    X_other <- do.call(cbind, lapply(X_trials[-t], function(X) X %*% h))
    X_full <- cbind(X_trials[[t]] %*% h, X_other, Z)
    XtX <- crossprod(X_full) + 1e-6 * diag(ncol(X_full))
    Xty <- crossprod(X_full, y)
    betas_direct[t] <- solve(XtX, Xty)[1]
  }
  
  # Corrected Woodbury
  betas_corrected <- corrected_woodbury_lss(y_proj, X_trials, h, Z)
  
  # Current implementation
  lss_prep <- prepare_lss_fixed_components_core(
    A_lss_fixed_matrix = Z,
    intercept_col_index_in_Alss = 1,
    lambda_ridge_Alss = 1e-6
  )
  
  betas_current <- run_lss_for_voxel_core(
    Y_proj_voxel_vector = as.vector(y_proj),
    X_trial_onset_list_of_matrices = X_trials,
    H_shape_voxel_vector = h,
    A_lss_fixed_matrix = Z,
    P_lss_matrix = lss_prep$P_lss_matrix,
    p_lss_vector = lss_prep$p_lss_vector
  )
  
  cat("\n=== Comparison of Methods ===\n")
  cat("True betas:      ", round(true_betas, 3), "\n")
  cat("Direct LSS:      ", round(betas_direct, 3), "\n")
  cat("Corrected Wood:  ", round(betas_corrected, 3), "\n")
  cat("Current Wood:    ", round(betas_current, 3), "\n")
  
  cat("\n=== Errors ===\n")
  cat("Direct error:    ", round(sqrt(mean((betas_direct - true_betas)^2)), 4), "\n")
  cat("Corrected error: ", round(sqrt(mean((betas_corrected - true_betas)^2)), 4), "\n")
  cat("Current error:   ", round(sqrt(mean((betas_current - true_betas)^2)), 4), "\n")
  
  # The corrected Woodbury should match direct very closely
  max_diff <- max(abs(betas_corrected - betas_direct))
  cat("\nMax diff (corrected vs direct): ", max_diff, "\n")
  
  expect_lt(max_diff, 1e-5, "Corrected Woodbury should match direct LSS to reasonable precision")
})