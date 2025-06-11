# Simplified LSS Implementation Using fmrilss Package
# All LSS functionality is routed through the external fmrilss implementation
# 
# This file provides the core LSS (Least Squares Separate) functionality for
# the manifold_hrf package. As of version 0.1.0, all LSS computations are
# delegated to the fmrilss package for improved performance and reliability.

#' Prepare LSS Fixed Components (Core)
#'
#' Pre-computes the pseudoinverse of fixed effect regressors for efficient
#' trial-wise LSS computation. This is used when projecting out confounds.
#'
#' @param A_fixed_regressors_matrix An n x q matrix of fixed regressors
#'   (e.g., intercept, motion parameters, drift terms)
#' @param intercept_column Index of intercept column (deprecated, ignored)
#' @param lambda_fixed Ridge penalty for fixed effects (scalar >= 0, default 0)
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{P_lss_matrix}: q x n pseudoinverse matrix for fixed effects
#'     \item \code{p_lss_vector}: length n projection vector for intercept (unused)
#'   }
#'
#' @export
prepare_lss_fixed_components_core <- function(A_fixed_regressors_matrix = NULL,
                                            intercept_column = NULL,
                                            lambda_fixed = 0,
                                            # Old argument names for backward compatibility
                                            A_lss_fixed_matrix = NULL,
                                            intercept_col_index_in_Alss = NULL,
                                            lambda_ridge_Alss = NULL) {
  # Handle old argument names
  if (!is.null(A_lss_fixed_matrix)) {
    A_fixed_regressors_matrix <- A_lss_fixed_matrix
  }
  if (!is.null(intercept_col_index_in_Alss)) {
    intercept_column <- intercept_col_index_in_Alss
  }
  if (!is.null(lambda_ridge_Alss)) {
    lambda_fixed <- lambda_ridge_Alss
  }
  
  # Handle old three-argument interface
  if (!is.null(intercept_column) && is.numeric(intercept_column) && 
      is.numeric(lambda_fixed) && lambda_fixed < 1) {
    # Old interface: second arg was intercept column, third was lambda
    lambda_fixed <- lambda_fixed
  } else if (is.numeric(intercept_column) && intercept_column >= 1) {
    # Second argument is actually lambda in new interface
    lambda_fixed <- intercept_column
  }
  if (!is.matrix(A_fixed_regressors_matrix)) {
    stop("A_fixed_regressors_matrix must be a matrix")
  }
  
  n <- nrow(A_fixed_regressors_matrix)
  q <- ncol(A_fixed_regressors_matrix)
  
  if (n < 2) {
    stop("A_fixed_regressors_matrix must have at least 2 rows (timepoints)")
  }
  
  if (q < 1) {
    stop("A_fixed_regressors_matrix must have at least 1 column")
  }
  
  if (q >= n) {
    stop(sprintf(
      "A_fixed_regressors_matrix has too many columns (%d) relative to rows (%d). Must have q < n.",
      q, n
    ))
  }
  
  if (!is.null(intercept_column)) {
    if (!is.numeric(intercept_column) || 
        length(intercept_column) != 1 ||
        intercept_column < 1 || 
        intercept_column > q ||
        intercept_column != round(intercept_column)) {
      stop(sprintf(
        "intercept_column must be an integer between 1 and %d", 
        q
      ))
    }
  }
  
  lambda_fixed <- .validate_and_standardize_lambda(lambda_fixed, "lambda_fixed")
  
  # Compute regularized pseudoinverse
  AtA <- crossprod(A_fixed_regressors_matrix)
  AtA_reg <- AtA + lambda_fixed * diag(q)
  
  # Check condition number and add jitter if needed
  eigvals <- eigen(AtA_reg, symmetric = TRUE, only.values = TRUE)$values
  min_eigval <- min(eigvals)
  max_eigval <- max(eigvals)
  condition_number <- max_eigval / max(min_eigval, .Machine$double.eps)
  
  if (condition_number > 1e8) {
    # Add small jitter to improve conditioning
    jitter <- 1e-6 * median(diag(AtA_reg))
    AtA_reg <- AtA_reg + jitter * diag(q)
    warning(sprintf(
      "Fixed regressor matrix is poorly conditioned (condition number = %.2e). Added jitter = %.2e to diagonal.",
      condition_number, jitter
    ))
  }
  
  # Use Cholesky decomposition for numerical stability
  P_lss <- tryCatch({
    AtA_reg_inv <- chol2inv(chol(AtA_reg))
    AtA_reg_inv %*% t(A_fixed_regressors_matrix)
  }, error = function(e) {
    warning("Cholesky decomposition failed, using SVD pseudoinverse")
    MASS::ginv(A_fixed_regressors_matrix)
  })
  
  # Step 5: Compute p_lss_vector for intercept handling
  if (!is.null(intercept_column)) {
    # Extract the row corresponding to the intercept
    p_lss_vector <- P_lss[intercept_column, , drop = TRUE]
  } else {
    # No intercept specified, use zero vector
    p_lss_vector <- rep(0, n)
  }
  
  list(
    P_lss_matrix = P_lss,
    p_lss_vector = p_lss_vector
  )
}

#' Reconstruct HRF Shapes from Manifold Coordinates (Core)
#'
#' Transforms smoothed manifold coordinates back to HRF shape space using
#' the reconstructor matrix from manifold construction.
#'
#' @param B_reconstructor_matrix The p x m reconstructor matrix from
#'   \code{get_manifold_basis_reconstructor_core}, where p is HRF length
#'   and m is manifold dimensionality
#' @param Xi_smoothed_matrix The m x V matrix of spatially smoothed manifold
#'   coordinates from Component 2, where V is number of voxels
#'
#' @return H_shapes_allvox_matrix A p x V matrix of reconstructed HRF shapes
#'
#' @export
reconstruct_hrf_shapes_core <- function(B_reconstructor_matrix,
                                      Xi_smoothed_matrix) {
  if (!is.matrix(B_reconstructor_matrix)) {
    stop("B_reconstructor_matrix must be a matrix")
  }
  
  if (!is.matrix(Xi_smoothed_matrix)) {
    stop("Xi_smoothed_matrix must be a matrix")
  }
  
  # Get dimensions
  p <- nrow(B_reconstructor_matrix)
  m_B <- ncol(B_reconstructor_matrix)
  m_Xi <- nrow(Xi_smoothed_matrix)
  V <- ncol(Xi_smoothed_matrix)
  
  # Check dimension compatibility
  if (m_B != m_Xi) {
    stop(sprintf(
      "Dimension mismatch: B_reconstructor_matrix has %d columns but Xi_smoothed_matrix has %d rows",
      m_B, m_Xi
    ))
  }
  
  # Reconstruct HRF shapes: H = B * Xi
  H_shapes_allvox_matrix <- B_reconstructor_matrix %*% Xi_smoothed_matrix
  
  # Ensure the result is a regular matrix (not Matrix class)
  if (inherits(H_shapes_allvox_matrix, "Matrix")) {
    H_shapes_allvox_matrix <- as.matrix(H_shapes_allvox_matrix)
  }
  
  return(H_shapes_allvox_matrix)
}

#' Run LSS for Single Voxel (Core)
#'
#' Computes trial-wise beta estimates using Least Squares Separate (LSS)
#' for a single voxel. This function delegates the actual LSS computation
#' to the fmrilss package for improved performance and consistency.
#' 
#' @details The LSS approach fits each trial separately, with all other trials
#' treated as nuisance regressors. This implementation uses fmrilss::lss()
#' internally. Note that the input data (Y_proj_voxel_vector) should already
#' be projected to remove confounds to avoid double projection.
#'
#' @param Y_proj_voxel_vector Projected data (n x 1) with confounds removed
#' @param X_trial_onset_list_of_matrices List of unprojected trial matrices (n x p each)
#' @param H_shape_voxel_vector HRF shape (p x 1)
#' @param A_lss_fixed_matrix Matrix of fixed regressors used during projection (n x q)
#' @param P_lss_matrix Precomputed matrix from \code{prepare_lss_fixed_components_core}
#' @param p_lss_vector Precomputed intercept projection vector (unused)
#' @return Vector of trial-wise betas
#' @export
run_lss_for_voxel_core <- function(Y_proj_voxel_vector,
                                  X_trial_onset_list_of_matrices,
                                  H_shape_voxel_vector,
                                  A_lss_fixed_matrix,
                                  P_lss_matrix,
                                  p_lss_vector) {
  
  n <- length(Y_proj_voxel_vector)
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  # Create convolved trial regressors
  C <- matrix(0, n, T_trials)
  for (t in seq_len(T_trials)) {
    C[, t] <- X_trial_onset_list_of_matrices[[t]] %*% H_shape_voxel_vector
  }
  
  # IMPORTANT: Y_proj_voxel_vector is already projected, so we should NOT 
  # pass confounds to fmrilss::lss to avoid double projection
  # Use fmrilss to compute LSS without confounds since data is pre-projected
  result <- fmrilss::lss(
    Y = matrix(Y_proj_voxel_vector, ncol = 1),
    X = C,
    Z = NULL,  # Don't pass confounds - data is already projected!
    method = "r_optimized"
  )
  
  as.vector(result)
}

#' Run LSS Across Voxels (Core)
#'
#' Computes trial-wise LSS estimates for all voxels using parallel processing.
#' Routes all computation through the fmrilss package.
#'
#' @param Y_proj_matrix n x V projected BOLD data matrix
#' @param X_trial_onset_list_of_matrices List of length T with n x p design matrices
#' @param H_shapes_allvox_matrix p x V matrix of HRF shapes
#' @param A_lss_fixed_matrix n x q matrix of fixed regressors
#' @param P_lss_matrix q x n matrix from \code{prepare_lss_fixed_components_core}
#' @param p_lss_vector Intercept projection vector (unused)
#' @param n_jobs Number of parallel workers for voxel-wise computation
#' @param ram_heuristic_GB_for_Rt RAM limit in gigabytes (currently unused)
#' 
#' @return A T x V matrix of trial-wise beta estimates
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
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  # Pre-compute all convolved trial regressors (n x T x V tensor)
  # For efficiency, we can compute C = X * H for all voxels at once
  all_C <- array(0, dim = c(n, T_trials, V))
  for (t in seq_len(T_trials)) {
    # X_t is n x p, H is p x V, result is n x V
    all_C[, t, ] <- X_trial_onset_list_of_matrices[[t]] %*% H_shapes_allvox_matrix
  }
  
  # Process voxels in parallel or serial
  voxel_fun <- function(v) {
    tryCatch({
      # Extract the trial regressors for this voxel
      C_v <- all_C[, , v]  # n x T matrix
      
      # Use fmrilss to compute LSS
      # Note: Y_proj_matrix is already projected, so we don't pass Z
      result <- fmrilss::lss(
        Y = Y_proj_matrix[, v, drop = FALSE],
        X = C_v,
        Z = NULL,  # Data is already projected
        method = "r_optimized"
      )
      
      as.vector(result)
    }, error = function(e) {
      warning(sprintf("LSS failed for voxel %d: %s", v, e$message))
      rep(NA_real_, T_trials)
    })
  }
  
  res_list <- .parallel_lapply(seq_len(V), voxel_fun, n_jobs)
  do.call(cbind, res_list)
}


# Simple wrapper for single voxel with HRF convolution
#' Run LSS for Single Voxel with Simple Interface
#'
#' @param y_voxel Data vector for single voxel (n x 1)
#' @param X_trial_list List of trial matrices (n x p each)
#' @param h_voxel HRF shape vector (p x 1)
#' @param TR Repetition time in seconds (currently unused)
#' @return List with beta_trials vector
#' @export
run_lss_for_voxel <- function(y_voxel, X_trial_list, h_voxel, TR = 2) {
  n <- length(y_voxel)
  T_trials <- length(X_trial_list)
  
  # Validate inputs
  if (!all(is.finite(y_voxel))) {
    warning("Non-finite values in y_voxel")
    return(list(beta_trials = rep(NA_real_, T_trials)))
  }
  
  if (!all(is.finite(h_voxel))) {
    warning("Non-finite values in h_voxel")
    return(list(beta_trials = rep(NA_real_, T_trials)))
  }
  
  # Create convolved regressors
  C <- matrix(0, n, T_trials)
  for (t in seq_len(T_trials)) {
    C[, t] <- X_trial_list[[t]] %*% h_voxel
  }
  
  # Check if C has valid variance
  col_vars <- apply(C, 2, var)
  if (any(col_vars < .Machine$double.eps)) {
    warning("Some trial regressors have zero variance - likely due to trial timing at data edges")
    # Still try to compute but results may be NA for zero-variance columns
  }
  
  # Use fmrilss with intercept only
  result <- tryCatch({
    fmrilss::lss(
      Y = matrix(y_voxel, ncol = 1),
      X = C,
      Z = matrix(1, n, 1),  # Intercept only
      method = "r_optimized"
    )
  }, error = function(e) {
    warning("LSS computation failed: ", e$message)
    rep(NA_real_, T_trials)
  })
  
  list(beta_trials = as.vector(result))
}



# Removed internal adapter functions - all LSS computation now goes through fmrilss::lss

#' Fast LSS implementation using fmrilss
#' @keywords internal
lss_fast <- function(Y, dmat_base, dmat_ran, dmat_fixed = NULL) {
  if (!is.matrix(Y)) stop("Y must be a matrix")
  
  # Combine fixed effects
  Z <- if (is.null(dmat_fixed)) dmat_base else cbind(dmat_base, dmat_fixed)
  
  # Use fmrilss to compute LSS betas
  result <- fmrilss::lss(
    Y = Y,
    X = dmat_ran,
    Z = Z,
    method = "r_optimized"
  )
  
  return(result)
}





# Validation functions
#' Validate Design Matrix List
#' @keywords internal
validate_design_matrix_list <- function(X_list, n_timepoints) {
  if (!is.list(X_list)) {
    stop("Design matrix must be a list")
  }
  k <- length(X_list)
  if (k == 0) {
    stop("Design matrix list cannot be empty")
  }
  p <- ncol(X_list[[1]])
  for (i in seq_len(k)) {
    if (!is.matrix(X_list[[i]])) {
      stop(sprintf("Design matrix %d is not a matrix", i))
    }
    if (nrow(X_list[[i]]) != n_timepoints) {
      stop(sprintf("Design matrix %d has wrong number of rows", i))
    }
    if (ncol(X_list[[i]]) != p) {
      stop(sprintf("Design matrix %d has inconsistent columns", i))
    }
  }
  list(k = k, p = p)
}

#' Validate HRF Shape Matrix
#' @keywords internal
validate_hrf_shape_matrix <- function(H_matrix, n_timepoints, n_voxels) {
  if (!is.matrix(H_matrix)) {
    stop("HRF shape matrix must be a matrix")
  }
  if (nrow(H_matrix) != n_timepoints) {
    stop("HRF shape matrix has wrong number of rows")
  }
  if (ncol(H_matrix) != n_voxels) {
    stop("HRF shape matrix has wrong number of columns")
  }
}

# User-facing wrapper functions
#' Run LSS Voxel Loop
#' @export
run_lss_voxel_loop <- function(Y_matrix, X_trial_list, H_matrix, 
                               Z_confounds = NULL, lambda = 0) {
  # Convert to core interface
  n <- nrow(Y_matrix)
  A_fixed <- if (is.null(Z_confounds)) matrix(1, n, 1) else Z_confounds
  
  # Prepare fixed components
  lss_prep <- prepare_lss_fixed_components_core(A_fixed, lambda)
  
  # Run core function
  run_lss_voxel_loop_core(
    Y_proj_matrix = Y_matrix,
    X_trial_onset_list_of_matrices = X_trial_list,
    H_shapes_allvox_matrix = H_matrix,
    A_lss_fixed_matrix = A_fixed,
    P_lss_matrix = lss_prep$P_lss_matrix,
    p_lss_vector = lss_prep$p_lss_vector,
    n_jobs = 1
  )
}

#' Validate and Standardize Lambda
#' @keywords internal
.validate_and_standardize_lambda <- function(lambda, param_name) {
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) {
    stop(param_name, " must be a non-negative scalar")
  }
  as.numeric(lambda)
}