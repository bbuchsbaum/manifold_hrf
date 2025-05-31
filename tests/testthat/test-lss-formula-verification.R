# Verify the LSS formula implementation against reference
library(testthat)

# Helper function: Solve for A^+ using Cholesky
cholSolve <- function(A, B) {
  chol_A <- chol(A)
  backsolve(chol_A, backsolve(chol_A, B, transpose = TRUE))
}

# Reference implementation of Woodbury residualize
woodbury_residualize <- function(C, A, lambda_ridge = 0) {
  n <- nrow(A)
  m <- ncol(A)
  
  if (m == 0) return(C)
  
  AtA <- crossprod(A)
  if (lambda_ridge != 0) {
    AtA <- AtA + lambda_ridge * diag(m)
  }
  
  # Solve for A^+ C where A^+ = (A'A + λI)^{-1}A'
  U <- cholSolve(AtA, crossprod(A, C))
  V <- C - A %*% U
  return(V)
}

test_that("Current LSS formula matches reference implementation", {
  
  set.seed(111)
  
  # Test setup
  n <- 100
  T_trials <- 6
  q <- 4  # confounds
  lambda <- 1e-6
  
  # Create test data
  C <- matrix(rnorm(n * T_trials), n, T_trials)  # Trial regressors
  A <- cbind(1, poly(1:n, degree = q-1))          # Confound matrix
  y <- rnorm(n) + C %*% rnorm(T_trials) + A %*% rnorm(q)
  
  # Project y
  P <- diag(n) - A %*% solve(crossprod(A) + lambda * diag(q)) %*% t(A)
  y_proj <- P %*% y
  
  # METHOD 1: Reference implementation
  V_ref <- woodbury_residualize(C, A, lambda)
  
  # Compute p_vec (projection of intercept)
  AtA_inv <- solve(crossprod(A) + lambda * diag(q))
  p_vec <- A %*% AtA_inv[, 1]  # Intercept projection uses all confounds
  
  # Reference formula
  pc_row <- drop(crossprod(p_vec, C))
  cv_row <- colSums(V_ref * V_ref)
  alpha_row <- numeric(T_trials)
  nz <- cv_row > 0
  alpha_row[nz] <- (1 - pc_row[nz]) / cv_row[nz]
  S_ref <- sweep(V_ref, 2, alpha_row, "*") + p_vec
  betas_ref <- drop(crossprod(S_ref, y_proj))
  
  # METHOD 2: Current implementation
  lss_prep <- prepare_lss_fixed_components_core(
    A_lss_fixed_matrix = A,
    intercept_col_index_in_Alss = 1,
    lambda_ridge_Alss = lambda
  )
  
  # Check if P_lss is correct
  # P_lss should be (A'A + λI)^{-1}A'
  P_lss_expected <- AtA_inv %*% t(A)
  cat("\nP_lss error: ", max(abs(lss_prep$P_lss_matrix - P_lss_expected)), "\n")
  
  # Manually compute V using current method
  U_v <- lss_prep$P_lss_matrix %*% C
  V_current <- C - A %*% U_v
  
  # Compare V matrices
  V_diff <- max(abs(V_current - V_ref))
  cat("V matrix difference: ", V_diff, "\n")
  expect_lt(V_diff, 1e-10, "V matrices should match")
  
  # Now run the full current implementation
  betas_current <- numeric(T_trials)
  for (t in 1:T_trials) {
    # Single trial at a time (as in run_lss_for_voxel_core)
    beta_t <- run_lss_for_voxel_core(
      Y_proj_voxel_vector = y_proj,
      X_trial_onset_list_of_matrices = lapply(1:T_trials, function(i) {
        if (i == t) matrix(C[, i], ncol = 1) else matrix(0, n, 1)
      }),
      H_shape_voxel_vector = 1,  # Simple case
      A_lss_fixed_matrix = A,
      P_lss_matrix = lss_prep$P_lss_matrix,
      p_lss_vector = lss_prep$p_lss_vector
    )
    betas_current[t] <- beta_t[t]
  }
  
  cat("\nBetas comparison:\n")
  cat("Reference:  ", round(betas_ref, 4), "\n")
  cat("Current:    ", round(betas_current, 4), "\n")
  cat("Difference: ", round(betas_current - betas_ref, 6), "\n")

  expect_equal(betas_current, betas_ref, tolerance = 1e-6)
})


test_that("LSS formula gives correct results for known problem", {
  
  # Create a simple problem where we know the answer
  set.seed(222)
  
  n <- 50
  T_trials <- 3
  
  # Simple design: non-overlapping trials
  C <- matrix(0, n, T_trials)
  C[1:10, 1] <- 1
  C[21:30, 2] <- 1  
  C[41:50, 3] <- 1
  
  # Simple confounds
  A <- cbind(1, (1:n)/n)
  
  # True parameters
  true_betas <- c(2, -1, 1.5)
  true_conf <- c(5, -3)
  
  # Generate data
  y <- C %*% true_betas + A %*% true_conf + rnorm(n, sd = 0.1)
  
  # Project y
  lambda <- 1e-6
  AtA_inv <- solve(crossprod(A) + lambda * diag(ncol(A)))
  P <- diag(n) - A %*% AtA_inv %*% t(A)
  y_proj <- P %*% y
  
  # Apply LSS formula
  V <- woodbury_residualize(C, A, lambda)
  p_vec <- A %*% AtA_inv[, 1]
  
  pc_row <- drop(crossprod(p_vec, C))
  cv_row <- colSums(V * V)
  alpha_row <- (1 - pc_row) / pmax(cv_row, .Machine$double.eps)
  S <- sweep(V, 2, alpha_row, "*") + p_vec
  betas_lss <- drop(crossprod(S, y_proj))
  
  cat("\n=== Simple Non-overlapping Trials ===\n")
  cat("True betas: ", true_betas, "\n")
  cat("LSS betas:  ", round(betas_lss, 3), "\n")
  cat("Error:      ", round(betas_lss - true_betas, 3), "\n")
  
  # For non-overlapping trials, LSS should recover the true values well
  expect_lt(sqrt(mean((betas_lss - true_betas)^2)), 0.2,
            "LSS should recover non-overlapping trials accurately")
})


test_that("p_lss_vector computation is correct", {
  
  # The p_lss_vector should be the projection of the intercept column
  
  set.seed(333)
  n <- 30
  q <- 3
  
  # Confound matrix with intercept first
  A <- cbind(1, matrix(rnorm(n * (q-1)), n, q-1))
  
  # Compute p_lss_vector using current implementation
  lss_prep <- prepare_lss_fixed_components_core(
    A_lss_fixed_matrix = A,
    intercept_col_index_in_Alss = 1,
    lambda_ridge_Alss = 1e-6
  )
  
  # Manual computation
  AtA_inv <- solve(crossprod(A) + 1e-6 * diag(q))
  
  # Method 1: P * e_1 where e_1 is the first standard basis vector
  P_lss <- AtA_inv %*% t(A)
  p_manual1 <- P_lss[1, ]  # First row of P_lss
  
  # Method 2: Using the formula from the code
  # ginv(A)[intercept_col, ] which is (A'A)^{-1}A'[intercept_col, ]
  p_manual2 <- A %*% AtA_inv[, 1]
  
  cat("\n=== p_lss_vector computation ===\n")
  cat("Current implementation norm: ", sqrt(sum(lss_prep$p_lss_vector^2)), "\n")
  cat("Manual method 1 norm:        ", sqrt(sum(p_manual1^2)), "\n")
  cat("Manual method 2 norm:        ", sqrt(sum(p_manual2^2)), "\n")
  
  # The current implementation uses P_lss[intercept_col, ]
  p_current_check <- lss_prep$P_lss_matrix[1, ]
  
  cat("\nDifference (current vs P_lss[1,]): ", 
      max(abs(lss_prep$p_lss_vector - p_current_check)), "\n")
  
  expect_equal(lss_prep$p_lss_vector, p_current_check, tolerance = 1e-10,
               "p_lss_vector should be the first row of P_lss")
})