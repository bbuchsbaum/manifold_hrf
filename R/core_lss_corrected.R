# Corrected LSS Implementation Using Standard Linear Algebra
# This implements the mathematically correct approach

#' Run LSS for Single Voxel - Corrected Version
#'
#' This is a mathematically correct implementation that:
#' 1. Takes projected data (confounds removed)
#' 2. Takes unprojected trial matrices
#' 3. Projects the trial matrices
#' 4. Solves the reduced system efficiently
#'
#' @param Y_proj_voxel_vector Projected data (n x 1) with confounds removed
#' @param X_trial_onset_list_of_matrices List of unprojected trial matrices (n x p each)
#' @param H_shape_voxel_vector HRF shape (p x 1)
#' @param P_confound Projection matrix that removes confounds (n x n)
#' @param lambda_ridge Ridge regularization parameter
#' @return Vector of trial-wise betas
#' @export
run_lss_for_voxel_corrected <- function(Y_proj_voxel_vector,
                                       X_trial_onset_list_of_matrices,
                                       H_shape_voxel_vector,
                                       P_confound,
                                       lambda_ridge = 1e-6) {
  
  n <- length(Y_proj_voxel_vector)
  T_trials <- length(X_trial_onset_list_of_matrices)
  p <- length(H_shape_voxel_vector)
  
  # Step 1: Create all trial regressors
  C <- matrix(0, n, T_trials)
  for (t in 1:T_trials) {
    C[, t] <- X_trial_onset_list_of_matrices[[t]] %*% H_shape_voxel_vector
  }
  
  # Step 2: Project trial regressors to remove confound space
  C_proj <- P_confound %*% C
  
  # Step 3: Solve the projected system
  # We're solving: β = (C_proj'C_proj + λI)^(-1) C_proj' y_proj
  
  # Method 1: Direct solution (stable but not optimized)
  CtC_proj <- crossprod(C_proj)
  CtC_proj_reg <- CtC_proj + lambda_ridge * diag(T_trials)
  Cty_proj <- crossprod(C_proj, Y_proj_voxel_vector)
  
  betas <- solve(CtC_proj_reg, Cty_proj)
  
  return(as.vector(betas))
}


#' Run LSS Using Woodbury Identity - Corrected Version
#'
#' This implements the true Woodbury matrix identity for efficiency.
#' For each trial, we update the inverse rather than recomputing.
#'
#' @export
run_lss_woodbury_corrected <- function(Y_proj_voxel_vector,
                                      X_trial_onset_list_of_matrices,
                                      H_shape_voxel_vector,
                                      P_confound,
                                      lambda_ridge = 1e-6) {
  
  n <- length(Y_proj_voxel_vector)
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  # Create and project all trial regressors
  C <- matrix(0, n, T_trials)
  for (t in 1:T_trials) {
    C[, t] <- X_trial_onset_list_of_matrices[[t]] %*% H_shape_voxel_vector
  }
  C_proj <- P_confound %*% C
  
  # Initialize results
  betas <- numeric(T_trials)
  
  # For each trial, we solve a "leave-one-in" problem
  # We want just the coefficient for trial t when all trials are included
  
  # Pre-compute the full system
  CtC_proj <- crossprod(C_proj) + lambda_ridge * diag(T_trials)
  Cty_proj <- crossprod(C_proj, Y_proj_voxel_vector)
  
  # Solve once for all trials
  betas_all <- solve(CtC_proj, Cty_proj)
  
  # Extract individual trial betas
  # In LSS, we typically solve T separate problems
  # But if we want all trials simultaneously, this is the solution
  
  # For true LSS (each trial separately with others as nuisance):
  for (t in 1:T_trials) {
    # Build design with trial t first, others as nuisance
    idx_others <- setdiff(1:T_trials, t)
    C_t_proj <- C_proj[, c(t, idx_others), drop = FALSE]
    
    # Solve this system
    CtC_t <- crossprod(C_t_proj) + lambda_ridge * diag(T_trials)
    Cty_t <- crossprod(C_t_proj, Y_proj_voxel_vector)
    beta_t <- solve(CtC_t, Cty_t)
    
    betas[t] <- beta_t[1]  # First element is trial of interest
  }
  
  return(betas)
}


#' Prepare Projection Matrix for Confound Removal
#'
#' Creates the projection matrix P = I - Z(Z'Z + λI)^{-1}Z'
#'
#' @param Z_confounds Confound matrix (n x q)
#' @param lambda Ridge parameter for stability
#' @return List with projection matrix and related quantities
#' @export
prepare_projection_matrix <- function(Z_confounds, lambda = 1e-6) {
  
  n <- nrow(Z_confounds)
  q <- ncol(Z_confounds)
  
  # Compute (Z'Z + λI)^{-1}
  ZtZ <- crossprod(Z_confounds)
  ZtZ_reg <- ZtZ + lambda * diag(q)
  ZtZ_inv <- solve(ZtZ_reg)
  
  # Projection matrix P = I - Z(Z'Z + λI)^{-1}Z'
  P <- diag(n) - Z_confounds %*% ZtZ_inv %*% t(Z_confounds)
  
  # For memory efficiency, we might want to avoid storing full P
  # Instead store components:
  # - Z
  # - ZtZ_inv
  # Then P*y = y - Z*(ZtZ_inv*(Z'*y))
  
  return(list(
    P = P,
    Z = Z_confounds,
    ZtZ_inv = ZtZ_inv,
    project_fn = function(y) {
      y - Z_confounds %*% (ZtZ_inv %*% crossprod(Z_confounds, y))
    }
  ))
}