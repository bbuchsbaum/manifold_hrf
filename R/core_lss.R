# Corrected LSS Implementation Using Standard Linear Algebra
# This implements the mathematically correct approach

#' Run LSS for Single Voxel - Simple Interface
#'
#' Simplified interface for trial-wise LSS that handles projection internally
#'
#' @param y_voxel Data vector for single voxel (n x 1)
#' @param X_trial_list List of trial matrices (n x p each)
#' @param h_voxel HRF shape vector (p x 1)
#' @param TR Repetition time in seconds
#' @param lambda Ridge regularization parameter
#' @return List with beta_trials vector
#' @export
run_lss_for_voxel_corrected <- function(y_voxel,
                                       X_trial_list,
                                       h_voxel,
                                       TR = 2,
                                       lambda = 1e-6) {
  
  n <- length(y_voxel)
  T_trials <- length(X_trial_list)
  p <- length(h_voxel)
  
  # Create convolved regressors
  C <- matrix(0, n, T_trials)
  for (t in 1:T_trials) {
    X_trial <- X_trial_list[[t]]
    # Convolve with HRF
    regressor <- X_trial %*% h_voxel
    C[, t] <- regressor
  }

  # Check for rank deficiency
  qr_C <- qr(C)
  if (qr_C$rank < T_trials) {
    warning(
      sprintf(
        "Trial regressors are rank deficient for this voxel: rank %d < %d",
        qr_C$rank,
        T_trials
      )
    )
  }
  
  # Ridge regression
  CTC <- t(C) %*% C
  diag(CTC) <- diag(CTC) + lambda
  beta_trials <- solve(CTC) %*% t(C) %*% y_voxel
  
  return(list(beta_trials = as.vector(beta_trials)))
}


#' Run LSS for Single Voxel - Full Interface
#'
#' Complete interface that takes projected data and projection matrices
#'
#' @param Y_proj_voxel_vector Projected data (n x 1) with confounds removed
#' @param X_trial_onset_list_of_matrices List of unprojected trial matrices (n x p each)
#' @param H_shape_voxel_vector HRF shape (p x 1)
#' @param P_confound Projection matrix that removes confounds (n x n)
#' @param lambda_ridge Ridge regularization parameter
#' @return Vector of trial-wise betas
#' @export
run_lss_for_voxel_corrected_full <- function(Y_proj_voxel_vector,
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

  # Check for rank deficiency
  qr_C <- qr(C)
  if (qr_C$rank < T_trials) {
    warning(
      sprintf(
        "Trial regressors are rank deficient for this voxel: rank %d < %d",
        qr_C$rank,
        T_trials
      )
    )
  }
  
  # Step 2: Project trial regressors to remove confound space
  C_proj <- P_confound %*% C
  
  # Step 3: Solve the projected system
  # We're solving: β = (C_proj'C_proj + λI)^(-1) C_proj' y_proj
  
  # Method 1: Direct solution (stable but not optimized)
  CtC_proj <- crossprod(C_proj)
  CtC_proj_reg <- CtC_proj + lambda_ridge * diag(T_trials)
  
  beta_trials <- solve(CtC_proj_reg) %*% crossprod(C_proj, Y_proj_voxel_vector)
  
  return(as.vector(beta_trials))
}


#' Run Woodbury LSS for Single Voxel
#'
#' Memory-efficient implementation using Woodbury matrix identity
#'
#' @param Y_proj_voxel_vector Projected data vector (n x 1)
#' @param X_trial_onset_list_of_matrices List of trial design matrices
#' @param H_shape_voxel_vector HRF shape vector (p x 1)
#' @param P_confound Projection matrix used to remove confound effects
#'   (n x n). If \code{NULL}, no projection is applied.
#' @param lambda_ridge Ridge regularization parameter
#' @return Vector of trial-wise betas
#' @export
run_lss_woodbury_corrected <- function(Y_proj_voxel_vector,
                                      X_trial_onset_list_of_matrices,
                                      H_shape_voxel_vector,
                                      P_confound = NULL,
                                      lambda_ridge = 1e-6) {
  
  n <- length(Y_proj_voxel_vector)
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  if (is.null(P_confound)) {
    # Identity projection if none supplied
    P_confound <- diag(n)
  }

  beta_trials <- numeric(T_trials)
  
  # For each trial, solve the LSS problem
  for (t in 1:T_trials) {
    
    # Create design matrix for this trial
    X_t <- X_trial_onset_list_of_matrices[[t]]
    C_t <- X_t %*% H_shape_voxel_vector
    C_t <- P_confound %*% C_t
    
    # Create design matrix for all other trials
    other_trials <- setdiff(1:T_trials, t)
    if (length(other_trials) > 0) {
      C_others <- matrix(0, n, length(other_trials))
      for (i in seq_along(other_trials)) {
        trial_idx <- other_trials[i]
        X_other <- X_trial_onset_list_of_matrices[[trial_idx]]
        C_others[, i] <- P_confound %*% (X_other %*% H_shape_voxel_vector)
      }
      
      # Full design matrix
      X_full <- cbind(C_t, C_others)
    } else {
      X_full <- matrix(C_t, ncol = 1)
    }
    
    # Solve regularized system
    XTX <- crossprod(X_full)
    diag(XTX) <- diag(XTX) + lambda_ridge
    
    beta_full <- solve(XTX) %*% crossprod(X_full, Y_proj_voxel_vector)
    beta_trials[t] <- beta_full[1]  # First coefficient is for trial t
  }
  
  return(beta_trials)
}


#' Prepare Projection Matrix
#'
#' Creates a projection matrix that removes confound space
#'
#' @param Z_confounds Confound matrix (n x q)
#' @param lambda Ridge regularization parameter
#' @return Projection matrix (n x n)
#' @export
prepare_projection_matrix <- function(Z_confounds, lambda = 1e-6) {
  
  n <- nrow(Z_confounds)
  q <- ncol(Z_confounds)
  
  if (q == 0) {
    # No confounds - return identity
    return(diag(n))
  }
  
  # Regularized projection: P = I - Z(Z'Z + λI)^(-1)Z'
  ZTZ <- crossprod(Z_confounds)
  ZTZ_reg <- ZTZ + lambda * diag(q)
  
  # Use Woodbury identity for stability
  ZTZ_inv <- solve(ZTZ_reg)
  P_confound <- diag(n) - Z_confounds %*% ZTZ_inv %*% t(Z_confounds)
  
  return(P_confound)
}
#' Run LSS Across Voxels (Core)
#'
#' Computes trial-wise LSS estimates for all voxels using the
#' Woodbury-based implementation.
#'
#' @param Y_proj_matrix n x V projected BOLD data matrix.
#' @param X_trial_onset_list_of_matrices List of length T with n x p design matrices.
#' @param H_shapes_allvox_matrix p x V matrix of HRF shapes.
#' @param A_lss_fixed_matrix n x q matrix of fixed regressors.
#' @param lambda_ridge Ridge regularization applied when projecting out confounds.
#' @param n_jobs Number of parallel workers for voxel-wise computation.
#' @param ram_heuristic_GB_for_Rt RAM limit in gigabytes for optional
#'   precomputation of trial regressors. Precomputation is used when the
#'   estimated memory footprint \code{T * V * 8 / 1e9} is below this limit.
#' 
#' @return A T x V matrix of trial-wise beta estimates.
#' @export
run_lss_voxel_loop_core <- function(Y_proj_matrix,
                                   X_trial_onset_list_of_matrices,
                                   H_shapes_allvox_matrix,
                                   A_lss_fixed_matrix,
                                   lambda_ridge = 1e-6,
                                   n_jobs = 1,
                                   ram_heuristic_GB_for_Rt = 1.0) {

  if (!is.matrix(Y_proj_matrix)) {
    stop("Y_proj_matrix must be a matrix")
  }

  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)

  design_info <- validate_design_matrix_list(
    X_trial_onset_list_of_matrices,
    n_timepoints = n
  )
  T_trials <- design_info$k

  validate_hrf_shape_matrix(
    H_shapes_allvox_matrix,
    n_timepoints = design_info$p,
    n_voxels = V
  )

  P_confound <- prepare_projection_matrix(A_lss_fixed_matrix, lambda_ridge)

  precompute_Rt <- check_ram_feasibility(T_trials, V, ram_heuristic_GB_for_Rt)

  if (precompute_Rt) {
    R_t_allvox_list <- vector("list", T_trials)
    for (t in seq_len(T_trials)) {
      R_t_allvox_list[[t]] <- X_trial_onset_list_of_matrices[[t]] %*% H_shapes_allvox_matrix
    }
  }

  voxel_fun <- function(v) {
##<<<<<<< codex/add-memory-based-precomputation-to-run_lss_voxel_loop_core
    Y_proj_voxel_vector <- Y_proj_matrix[, v]

    if (precompute_Rt) {
      C_v <- matrix(0, n, T_trials)
      for (t in seq_len(T_trials)) {
        C_v[, t] <- R_t_allvox_list[[t]][, v]
      }

      qr_C <- qr(C_v)
      if (qr_C$rank < T_trials) {
        warning(
          sprintf(
            "Trial regressors are rank deficient for voxel %d: rank %d < %d",
            v, qr_C$rank, T_trials
          )
        )
      }

      C_proj <- P_confound %*% C_v
      CtC_proj <- crossprod(C_proj)
      CtC_proj_reg <- CtC_proj + lambda_ridge * diag(T_trials)
      beta_trials <- solve(CtC_proj_reg) %*% crossprod(C_proj, Y_proj_voxel_vector)
      as.vector(beta_trials)
    } else {
      run_lss_for_voxel_corrected_full(
        Y_proj_voxel_vector = Y_proj_matrix[, v],
        X_trial_onset_list_of_matrices = X_trial_onset_list_of_matrices,
        H_shape_voxel_vector = H_shapes_allvox_matrix[, v],
        P_confound = P_confound,
        lambda_ridge = lambda_ridge
      )
    }
##=======
    run_lss_woodbury_corrected(
      Y_proj_voxel_vector = Y_proj_matrix[, v],
      X_trial_onset_list_of_matrices = X_trial_onset_list_of_matrices,
      H_shape_voxel_vector = H_shapes_allvox_matrix[, v],
      P_confound = P_confound,
      lambda_ridge = lambda_ridge
    )
##>>>>>>> main
  }

  res_list <- .parallel_lapply(seq_len(V), voxel_fun, n_jobs)
  Beta_trial_allvox_matrix <- do.call(cbind, res_list)

  return(Beta_trial_allvox_matrix)
}
