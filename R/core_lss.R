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

  # Create convolved regressors
  C <- matrix(0, n, T_trials)
  for (t in seq_len(T_trials)) {
    C[, t] <- X_trial_list[[t]] %*% h_voxel
  }

  # Project out intercept
  X_base <- matrix(1, n, 1)
  P_base <- MASS::ginv(X_base)
  Q_base <- diag(n) - X_base %*% P_base
  y_res <- Q_base %*% y_voxel
  Q_dmat <- Q_base %*% C

  betas <- lss_compute_r(Q_dmat, as.matrix(y_res))
  list(beta_trials = as.vector(betas))
}


#' Compute LSS Betas for a Single Voxel
#'
#' Internal helper used by LSS functions to compute trial-wise betas using
#' the Woodbury formulation. Inputs should already be convolved with the
#' voxel-specific HRF.
#'
#' @param C_v Matrix of convolved trial regressors (n x T)
#' @param Y_v Numeric vector of length n with projected voxel data
#' @param A_fixed Fixed regressors matrix used during projection (n x q)
#' @param P_lss Precomputed matrix from `prepare_lss_fixed_components_core`
#'              with dimensions q x n
#' @param p_lss Precomputed intercept projection vector of length n
#' @return Numeric vector of length T containing LSS beta estimates
#' @keywords internal
.compute_lss_betas <- function(C_v, Y_v, A_fixed, P_lss, p_lss) {
  n <- nrow(C_v)
  P_confound <- diag(n) - A_fixed %*% P_lss
  Q_dmat <- P_confound %*% C_v
  lss_compute_r(Q_dmat, as.matrix(Y_v))
}


#' Run LSS for Single Voxel - Full Interface
#'
#' Complete interface that uses precomputed LSS components
#'
#' @param Y_proj_voxel_vector Projected data (n x 1) with confounds removed
#' @param X_trial_onset_list_of_matrices List of unprojected trial matrices (n x p each)
#' @param H_shape_voxel_vector HRF shape (p x 1)
#' @param A_lss_fixed_matrix Matrix of fixed regressors used during projection (n x q)
#' @param P_lss_matrix Precomputed matrix from \code{prepare_lss_fixed_components_core}
#'   with dimensions q x n
#' @param p_lss_vector Precomputed intercept projection vector of length n
#' @return Vector of trial-wise betas
#' @export
run_lss_for_voxel_corrected_full <- function(Y_proj_voxel_vector,
                                            X_trial_onset_list_of_matrices,
                                            H_shape_voxel_vector,
                                            A_lss_fixed_matrix,
                                            P_lss_matrix,
                                            p_lss_vector) {
  
  n <- length(Y_proj_voxel_vector)
  T_trials <- length(X_trial_onset_list_of_matrices)
  p <- length(H_shape_voxel_vector)
  
  # Step 1: Create all trial regressors
  C <- matrix(0, n, T_trials)
  for (t in seq_len(T_trials)) {
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
  
  .compute_lss_betas(
    C_v = C,
    Y_v = Y_proj_voxel_vector,
    A_fixed = A_lss_fixed_matrix,
    P_lss = P_lss_matrix,
    p_lss = p_lss_vector
  )
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

  # Build convolved trial regressors and project
  C <- matrix(0, n, T_trials)
  for (t in seq_len(T_trials)) {
    C[, t] <- X_trial_onset_list_of_matrices[[t]] %*% H_shape_voxel_vector
  }
  Q_dmat <- P_confound %*% C

  # Compute betas via numerically stable formula
  as.vector(lss_compute_r(Q_dmat, as.matrix(Y_proj_voxel_vector)))
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
#' @param P_lss_matrix q x n matrix from \code{prepare_lss_fixed_components_core}.
#' @param p_lss_vector Intercept projection vector from \code{prepare_lss_fixed_components_core}.
#' @param n_jobs Number of parallel workers for voxel-wise computation.
#' @param ram_heuristic_GB_for_Rt RAM limit in gigabytes (currently unused).
#' 
#' @return A T x V matrix of trial-wise beta estimates.
#' @export
run_lss_voxel_loop_core <- function(Y_proj_matrix,
                                   X_trial_onset_list_of_matrices,
                                   H_shapes_allvox_matrix,
                                   A_lss_fixed_matrix,
                                   P_lss_matrix,
                                   p_lss_vector,
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

  voxel_fun <- function(v) {
    run_lss_for_voxel_corrected_full(
      Y_proj_voxel_vector = Y_proj_matrix[, v],
      X_trial_onset_list_of_matrices = X_trial_onset_list_of_matrices,
      H_shape_voxel_vector = H_shapes_allvox_matrix[, v],
      A_lss_fixed_matrix = A_lss_fixed_matrix,
      P_lss_matrix = P_lss_matrix,
      p_lss_vector = p_lss_vector
    )
  }

  res_list <- .parallel_lapply(seq_len(V), voxel_fun, n_jobs)
  do.call(cbind, res_list)
}
