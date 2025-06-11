# Core Voxel-wise Manifold Fit Functions (Component 1)
# Implementation of MHRF-CORE-VOXFIT-01 through MHRF-CORE-VOXFIT-05

# Internal helper function for vectorized canonical correlation identifiability
.apply_identifiability_vectorized <- function(Xi_raw_matrix, Beta_raw_matrix,
                                             B_reconstructor_matrix, h_ref_shape_vector,
                                             scale_method, correlation_threshold,
                                             zero_tol, consistency_check) {
  
  m <- nrow(Xi_raw_matrix)
  V <- ncol(Xi_raw_matrix)
  k <- nrow(Beta_raw_matrix)
  p <- nrow(B_reconstructor_matrix)
  
  SCALE_OVERFLOW_THRESHOLD <- 1 / sqrt(.Machine$double.eps)
  
  # Step 1: Reconstruct all HRFs at once (p x V matrix)
  H_matrix <- B_reconstructor_matrix %*% Xi_raw_matrix
  
  # Step 2: Vectorized correlation with h_ref
  # Center h_ref
  h_ref_c <- h_ref_shape_vector - mean(h_ref_shape_vector)
  h_ref_norm <- sqrt(sum(h_ref_c^2))
  
  # For each column of H, compute correlation with h_ref
  # Using efficient column-wise operations
  H_means <- colMeans(H_matrix)
  H_centered <- sweep(H_matrix, 2, H_means, "-")
  H_norms <- sqrt(colSums(H_centered^2))
  
  # Compute correlations
  correlations <- as.vector(crossprod(h_ref_c, H_centered)) / (h_ref_norm * H_norms)
  correlations[is.na(correlations) | is.infinite(correlations)] <- 0
  
  # Step 3: Determine sign for all voxels
  sgn_vec <- sign(correlations)
  sgn_vec[sgn_vec == 0] <- 1
  
  # Fallback logic for low correlation (RMS rule)
  low_corr_idx <- abs(correlations) < correlation_threshold
  if (any(low_corr_idx)) {
    h_sums <- colSums(H_matrix[, low_corr_idx, drop = FALSE])
    fallback_sgn <- sign(h_sums)
    fallback_sgn[fallback_sgn == 0] <- 1
    sgn_vec[low_corr_idx] <- fallback_sgn
  }
  
  # Step 4: Apply sign to Xi and Beta
  Xi_signed <- sweep(Xi_raw_matrix, 2, sgn_vec, "*")
  Beta_signed <- sweep(Beta_raw_matrix, 2, sgn_vec, "*")
  
  # Step 5: Reconstruct HRFs with correct sign for scaling
  H_signed <- B_reconstructor_matrix %*% Xi_signed
  
  # Step 6: Apply scale normalization
  if (scale_method == "l2_norm") {
    # L2 norms of each column
    l2_norms <- sqrt(colSums(H_signed^2))
    # Handle zero norms
    zero_idx <- l2_norms < zero_tol
    scl_vec <- rep(1, V)
    scl_vec[!zero_idx] <- 1 / pmax(l2_norms[!zero_idx], .Machine$double.eps)
    
  } else if (scale_method == "max_abs_val") {
    # Maximum absolute value of each column
    max_abs_vals <- apply(abs(H_signed), 2, max)
    # Handle zero max values
    zero_idx <- max_abs_vals < zero_tol
    scl_vec <- rep(1, V)
    scl_vec[!zero_idx] <- 1 / pmax(max_abs_vals[!zero_idx], .Machine$double.eps)
    
  } else {  # "none"
    scl_vec <- rep(1, V)
    zero_idx <- rep(FALSE, V)
  }
  
  # Apply scaling
  Xi_out <- sweep(Xi_signed, 2, scl_vec, "*")
  Beta_out <- sweep(Beta_signed, 2, 1/scl_vec, "*")
  
  # Handle numerical overflow
  overflow_idx <- scl_vec > SCALE_OVERFLOW_THRESHOLD
  if (any(overflow_idx)) {
    Beta_out[, overflow_idx] <- 0
  }
  
  # Zero out voxels with no signal
  if (any(zero_idx)) {
    Xi_out[, zero_idx] <- 0
    Beta_out[, zero_idx] <- 0
  }
  
  # Step 7: Optional consistency check
  if (consistency_check) {
    # Reconstruct final HRFs
    H_final <- B_reconstructor_matrix %*% Xi_out
    
    # Recompute correlations with h_ref for final check
    H_final_means <- colMeans(H_final)
    H_final_centered <- sweep(H_final, 2, H_final_means, "-")
    H_final_norms <- sqrt(colSums(H_final_centered^2))
    
    final_correlations <- as.vector(crossprod(h_ref_c, H_final_centered)) / 
                         (h_ref_norm * H_final_norms)
    final_correlations[is.na(final_correlations)] <- 0
    
    # Flip voxels with negative correlation
    flip_idx <- final_correlations < 0
    if (any(flip_idx)) {
      Xi_out[, flip_idx] <- -Xi_out[, flip_idx]
      Beta_out[, flip_idx] <- -Beta_out[, flip_idx]
    }
  }
  
  list(
    Xi_ident_matrix = Xi_out,
    Beta_ident_matrix = Beta_out
  )
}

# Internal helper function for voxel processing in identifiability
.process_voxel_identifiability <- function(vx, Xi_raw, Beta_raw, B_reconstructor, 
                                         h_ref, config, Y_proj = NULL, 
                                         X_list = NULL, verbose = FALSE) {
  
  xi_vx <- Xi_raw[, vx]
  beta_vx <- Beta_raw[, vx]
  m <- length(xi_vx)
  k <- length(beta_vx)
  
  # Skip if no signal
  if (all(abs(xi_vx) < .Machine$double.eps)) {
    return(list(xi = rep(0, m), beta = rep(0, k)))
  }
  
  # Constants
  CORRELATION_THRESHOLD <- 1e-3
  SCALE_OVERFLOW_THRESHOLD <- 1 / sqrt(.Machine$double.eps)
  
  # Step 1: Determine sign
  sgn <- 1
  method_used <- config$sign_method
  
  if (config$sign_method == "canonical_correlation") {
    h_tmp <- B_reconstructor %*% xi_vx
    # Ensure both are vectors for correlation
    corr_ref <- tryCatch(cor(as.vector(h_tmp), as.vector(h_ref)),
                        warning = function(w) NA,
                        error = function(e) NA)
    if (is.na(corr_ref)) corr_ref <- 0

    sgn <- sign(corr_ref)
    if (sgn == 0) sgn <- 1

    # Fallback directly to RMS rule if correlation is too low
    if (abs(corr_ref) < CORRELATION_THRESHOLD) {
      sgn <- sign(sum(h_tmp))
      if (sgn == 0) sgn <- 1
      method_used <- "rms_fallback"
    }
  } else if (config$sign_method == "data_fit_correlation") {
    # Data fit method
    best_r2 <- -Inf
    best_sgn <- 1
    
    for (sg in c(1, -1)) {
      xi_tmp <- xi_vx * sg
      beta_tmp <- beta_vx * sg
      h_tmp <- B_reconstructor %*% xi_tmp
      
      X_design <- matrix(0, nrow(Y_proj), length(X_list))
      for (c in seq_along(X_list)) {
        X_design[, c] <- X_list[[c]] %*% h_tmp
      }
      
      y_pred <- X_design %*% beta_tmp
      y_true <- Y_proj[, vx]
      # Ensure both are vectors for correlation
      r2 <- tryCatch(cor(as.vector(y_pred), as.vector(y_true))^2, 
                    warning = function(w) NA, 
                    error = function(e) NA)
      
      if (!is.na(r2) && r2 > best_r2) {
        best_r2 <- r2
        best_sgn <- sg
      }
    }
    
    sgn <- best_sgn
    
    # Fallback if R² is too low
    if (best_r2 <= 0) {
      h_tmp <- B_reconstructor %*% xi_vx
      # Ensure both are vectors for correlation
      corr_ref <- tryCatch(cor(as.vector(h_tmp), as.vector(h_ref)), 
                          warning = function(w) NA, 
                          error = function(e) NA)
      
      if (!is.na(corr_ref) && abs(corr_ref) >= CORRELATION_THRESHOLD) {
        sgn <- sign(corr_ref)
        if (sgn == 0) sgn <- 1
        method_used <- "canonical_fallback"
      } else {
        sgn <- sign(sum(h_tmp))
        if (sgn == 0) sgn <- 1
        method_used <- "rms_fallback"
      }
    }
  }
  
  if (verbose) {
    message(sprintf("Voxel %d: sign method = %s", vx, method_used))
  }
  
  # Apply sign
  xi_signed <- xi_vx * sgn
  beta_signed <- beta_vx * sgn
  
  # Step 2: Reconstruct HRF for scaling
  hrf <- B_reconstructor %*% xi_signed
  
  # Step 3: Apply scale normalization
  if (config$scale_method == "l2_norm") {
    l2_norm <- sqrt(sum(hrf^2))
    if (l2_norm < config$zero_tol) {
      return(list(xi = rep(0, m), beta = rep(0, k)))
    }
    scl <- 1 / max(l2_norm, .Machine$double.eps)
  } else if (config$scale_method == "max_abs_val") {
    max_abs <- max(abs(hrf))
    if (max_abs < config$zero_tol) {
      return(list(xi = rep(0, m), beta = rep(0, k)))
    }
    scl <- 1 / max(max_abs, .Machine$double.eps)
  } else {  # "none"
    scl <- 1
  }
  
  # Apply scaling
  xi_out <- xi_signed * scl
  beta_out <- beta_signed / scl
  
  # Handle numerical overflow
  if (scl > SCALE_OVERFLOW_THRESHOLD) {
    beta_out <- rep(0, k)
  }
  
  # Step 4: Optional consistency check
  if (config$consistency_check) {
    hrf_check <- B_reconstructor %*% xi_out
    # Ensure both are vectors for correlation
    corr_check <- tryCatch(cor(as.vector(hrf_check), as.vector(h_ref)), 
                          warning = function(w) NA, 
                          error = function(e) NA)
    if (!is.na(corr_check) && corr_check < 0) {
      xi_out <- -xi_out
      beta_out <- -beta_out
    }
  }
  
  list(xi = xi_out, beta = beta_out)
}

#' Project Out Confounds (Core)
#'
#' Projects out confound variables from both the data matrix and design matrices
#' using QR decomposition.
#'
#' @param Y_data_matrix An n x V matrix of BOLD data (n timepoints, V voxels)
#' @param X_list_of_matrices A list of k matrices, each n x p (design matrices 
#'   for k conditions)
#' @param Z_confounds_matrix An n x q_confound matrix of confound regressors, 
#'   or NULL if no confounds
#'   
#' @return A list containing:
#'   \itemize{
#'     \item \code{Y_proj_matrix}: The n x V projected data matrix
#'     \item \code{X_list_proj_matrices}: List of k projected design matrices
#'   }
#'   
#' @details This function implements Component 1, Step 1 of the M-HRF-LSS pipeline.
#'   It uses SVD decomposition to robustly detect rank and project out confound 
#'   regressors from both the data and design matrices. If \code{Z_confounds_matrix} 
#'   is \code{NULL}, the original matrices are returned unchanged. Rank-deficient 
#'   confound matrices are automatically reduced to their independent columns with a
#'   warning. The SVD approach provides superior numerical stability compared to
#'   QR decomposition when dealing with near-collinear confounds. Missing values 
#'   are not allowed.
#'   
#' @examples
#' \dontrun{
#' # Create example data
#' n <- 200  # timepoints
#' V <- 100  # voxels
#' p <- 30   # HRF length
#' k <- 3    # conditions
#' 
#' Y_data <- matrix(rnorm(n * V), n, V)
#' X_list <- lapply(1:k, function(i) matrix(rnorm(n * p), n, p))
#' Z_confounds <- cbind(1, poly(1:n, degree = 3))  # intercept + polynomial trends
#' 
#' # Project out confounds
#' result <- project_out_confounds_core(Y_data, X_list, Z_confounds)
#' }
#' 
#' @export
project_out_confounds_core <- function(Y_data_matrix,
                                      X_list_of_matrices,
                                      Z_confounds_matrix = NULL) {
  
  # Input validation
  if (!is.matrix(Y_data_matrix)) {
    stop("Y_data_matrix must be a matrix")
  }

  if (anyNA(Y_data_matrix)) {
    stop("Y_data_matrix must not contain NA values")
  }

  n <- nrow(Y_data_matrix)

  validate_design_matrix_list(X_list_of_matrices, n_timepoints = n)
  
  # If no confounds, return original matrices
  if (is.null(Z_confounds_matrix)) {
    return(list(
      Y_proj_matrix = Y_data_matrix,
      X_list_proj_matrices = X_list_of_matrices
    ))
  }
  
  # Validate confounds matrix
  validate_confounds_matrix(Z_confounds_matrix, n_timepoints = n)
  
  # Check for rank deficiency in confounds
  if (ncol(Z_confounds_matrix) >= n) {
    stop("Z_confounds_matrix has too many columns (must be less than number of timepoints)")
  }

  # Step 1: Compute orthonormal basis of confounds (via SVD for robust rank detection)
  svd_Z <- svd(Z_confounds_matrix)
  tol_svd <- max(dim(Z_confounds_matrix)) * max(svd_Z$d) * .Machine$double.eps
  rank_Z <- sum(svd_Z$d > tol_svd)
  if (rank_Z < ncol(Z_confounds_matrix)) {
    warning("Z_confounds_matrix is rank deficient; using independent columns only")
  }
  Q_Z <- svd_Z$u[, seq_len(rank_Z), drop = FALSE]
  
  # Step 2: Project out confounds from Y
  # Y_proj = Y - Q_Z * Q_Z' * Y
  Y_proj_matrix <- Y_data_matrix - Q_Z %*% crossprod(Q_Z, Y_data_matrix)
  
  # Step 3: Project out confounds from each X matrix
  X_list_proj_matrices <- lapply(X_list_of_matrices, function(X) {
    X - Q_Z %*% crossprod(Q_Z, X)
  })
  
  # Return projected matrices
  return(list(
    Y_proj_matrix = Y_proj_matrix,
    X_list_proj_matrices = X_list_proj_matrices
  ))
}

#' Transform Designs to Manifold Basis (Core)
#'
#' Transforms design matrices from HRF space to manifold coordinate space.
#'
#' @param X_condition_list_proj_matrices List of k projected design matrices 
#'   (n x p each), where n is timepoints and p is HRF length
#' @param B_reconstructor_matrix The p x m HRF reconstructor matrix from 
#'   get_manifold_basis_reconstructor_core
#'   
#' @return Z_list_of_matrices List of k matrices (n x m each), where m is the 
#'   manifold dimensionality
#'   
#' @details This function implements Component 1, Step 2 of the M-HRF-LSS pipeline.
#'   It transforms each condition's design matrix from the original HRF space
#'   (p dimensions) to the lower-dimensional manifold space (m dimensions) by
#'   matrix multiplication with the reconstructor basis.
#'   
#' @examples
#' \dontrun{
#' # Create example data
#' n <- 200  # timepoints
#' p <- 30   # HRF length
#' m <- 5    # manifold dimensions
#' k <- 3    # conditions
#' 
#' # Example design matrices and reconstructor
#' X_list <- lapply(1:k, function(i) matrix(rnorm(n * p), n, p))
#' B_reconstructor <- matrix(rnorm(p * m), p, m)
#' 
#' # Transform to manifold basis
#' Z_list <- transform_designs_to_manifold_basis_core(X_list, B_reconstructor)
#' # Each Z_list[[i]] is now n x m instead of n x p
#' }
#' 
#' @export
transform_designs_to_manifold_basis_core <- function(X_condition_list_proj_matrices,
                                                    B_reconstructor_matrix) {
  
  # Input validation
  if (!is.list(X_condition_list_proj_matrices)) {
    stop("X_condition_list_proj_matrices must be a list")
  }
  
  if (!is.matrix(B_reconstructor_matrix)) {
    stop("B_reconstructor_matrix must be a matrix")
  }
  
  if (length(X_condition_list_proj_matrices) == 0) {
    stop("X_condition_list_proj_matrices cannot be empty")
  }
  
  # Get dimensions
  p <- nrow(B_reconstructor_matrix)
  m <- ncol(B_reconstructor_matrix)
  
  # Check that all X matrices have compatible dimensions
  for (i in seq_along(X_condition_list_proj_matrices)) {
    if (!is.matrix(X_condition_list_proj_matrices[[i]])) {
      stop(sprintf("X_condition_list_proj_matrices[[%d]] must be a matrix", i))
    }
    
    if (ncol(X_condition_list_proj_matrices[[i]]) != p) {
      stop(sprintf(
        "X_condition_list_proj_matrices[[%d]] has %d columns but B_reconstructor_matrix has %d rows",
        i, ncol(X_condition_list_proj_matrices[[i]]), p
      ))
    }
  }
  
  # Transform each design matrix to manifold basis
  # Z_i = X_i %*% B_reconstructor_matrix
  Z_list_of_matrices <- lapply(X_condition_list_proj_matrices, function(X_i) {
    X_i %*% B_reconstructor_matrix
  })
  
  return(Z_list_of_matrices)
}

#' Solve GLM for Gamma Coefficients (Core)
#'
#' Solves the GLM to estimate gamma coefficients for all voxels simultaneously.
#'
#' @param Z_list_of_matrices List of k design matrices in manifold basis (n x m each),
#'   where n is timepoints, m is manifold dimensions, and k is number of conditions
#' @param Y_proj_matrix The n x V projected data matrix, where V is number of voxels
#' @param lambda_gamma Ridge penalty parameter (scalar, typically small like 0.01)
#' @param orthogonal_approx_flag Boolean for orthogonal approximation. If TRUE,
#'   zeros out off-diagonal blocks in the design matrix cross-product, treating
#'   conditions as approximately orthogonal.
#'   
#' @return Gamma_coeffs_matrix A (km) x V matrix of gamma coefficients, where
#'   rows are organized as condition1_dim1, ..., condition1_dimM, 
#'   condition2_dim1, ..., condition2_dimM, etc.
#'   
#' @details This function implements Component 1, Step 3 of the M-HRF-LSS pipeline.
#'   It combines all condition design matrices into a single large design matrix
#'   and solves the ridge regression problem for all voxels at once. The optional
#'   orthogonal approximation can improve computational efficiency and stability
#'   when conditions are expected to be relatively independent.
#'   
#' @examples
#' \dontrun{
#' # Create example data
#' n <- 200  # timepoints
#' m <- 5    # manifold dimensions
#' k <- 3    # conditions
#' V <- 100  # voxels
#' 
#' # Design matrices in manifold basis
#' Z_list <- lapply(1:k, function(i) matrix(rnorm(n * m), n, m))
#' Y_proj <- matrix(rnorm(n * V), n, V)
#' 
#' # Solve GLM
#' gamma <- solve_glm_for_gamma_core(Z_list, Y_proj, lambda_gamma = 0.01)
#' # gamma is (k*m) x V
#' }
#' 
#' @export
solve_glm_for_gamma_core <- function(Z_list_of_matrices,
                                    Y_proj_matrix,
                                    lambda_gamma,
                                    orthogonal_approx_flag = FALSE) {
  
  # Input validation
  if (!is.list(Z_list_of_matrices)) {
    stop("Z_list_of_matrices must be a list")
  }
  
  if (length(Z_list_of_matrices) == 0) {
    stop("Z_list_of_matrices cannot be empty")
  }
  
  if (!is.matrix(Y_proj_matrix)) {
    stop("Y_proj_matrix must be a matrix")
  }
  
  lambda_gamma <- .validate_and_standardize_lambda(lambda_gamma, "lambda_gamma")
  
  if (!is.logical(orthogonal_approx_flag) || length(orthogonal_approx_flag) != 1) {
    stop("orthogonal_approx_flag must be a single logical value")
  }
  
  # Get dimensions
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  k <- length(Z_list_of_matrices)
  
  # Check first Z matrix to get m
  if (!is.matrix(Z_list_of_matrices[[1]])) {
    stop("Z_list_of_matrices[[1]] must be a matrix")
  }
  m <- ncol(Z_list_of_matrices[[1]])
  
  # Validate all Z matrices
  for (i in seq_along(Z_list_of_matrices)) {
    if (!is.matrix(Z_list_of_matrices[[i]])) {
      stop(sprintf("Z_list_of_matrices[[%d]] must be a matrix", i))
    }
    if (nrow(Z_list_of_matrices[[i]]) != n) {
      stop(sprintf("Z_list_of_matrices[[%d]] must have %d rows to match Y_proj_matrix", i, n))
    }
    if (ncol(Z_list_of_matrices[[i]]) != m) {
      stop(sprintf("All Z matrices must have the same number of columns. Z[[%d]] has %d columns but Z[[1]] has %d", 
                  i, ncol(Z_list_of_matrices[[i]]), m))
    }
  }
  
  # Step 1: Form combined design matrix X_tilde by column-binding all Z matrices
  # X_tilde is n x (k*m)
  X_tilde <- do.call(cbind, Z_list_of_matrices)
  
  # Step 2: Compute X'X
  XtX <- crossprod(X_tilde)  # (k*m) x (k*m)
  
  # Step 3: Apply orthogonal approximation if requested
  if (orthogonal_approx_flag) {
    # Zero out off-diagonal m x m blocks
    # The matrix is organized as k x k blocks, each of size m x m
    # We want to keep only the diagonal blocks
    
    # Create a block diagonal mask
    XtX_approx <- matrix(0, nrow = k * m, ncol = k * m)
    
    for (i in 1:k) {
      # Indices for block i
      idx <- ((i - 1) * m + 1):(i * m)
      # Copy the diagonal block
      XtX_approx[idx, idx] <- XtX[idx, idx]
    }
    
    XtX <- XtX_approx
  }
  
  # Step 4: Add ridge penalty
  XtX_reg <- XtX + lambda_gamma * diag(k * m)
  
  # Check condition number of regularized matrix
  cond_num <- kappa(XtX_reg, exact = FALSE)
  if (cond_num > 1e10) {
    warning(sprintf("XtX_reg has high condition number (%.2e). Consider increasing lambda_gamma.", cond_num))
  }
  
  # Step 5: Compute X'Y
  XtY <- crossprod(X_tilde, Y_proj_matrix)  # (k*m) x V
  
  # Step 6: Solve for gamma coefficients
  # Gamma = (X'X + lambda*I)^(-1) * X'Y
  # Use qr.solve for better numerical stability
  Gamma_coeffs_matrix <- qr.solve(XtX_reg, XtY)
  
  return(Gamma_coeffs_matrix)
}

#' Extract Xi and Beta via SVD (Core)
#'
#' Extracts raw manifold coordinates (Xi) and condition amplitudes (Beta) from
#' gamma coefficients using singular value decomposition.
#'
#' @param Gamma_coeffs_matrix The (km) x V coefficient matrix from 
#'   solve_glm_for_gamma_core, where rows are organized by condition and 
#'   manifold dimension
#' @param m_manifold_dim Manifold dimensionality (m)
#' @param k_conditions Number of conditions (k)
#' @param n_jobs Number of parallel jobs for voxel processing (default 1).
#' @param log_fun Optional logging function to report quality metrics. If NULL,
#'   no logging is performed.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{Xi_raw_matrix}: m x V matrix of raw manifold coordinates
#'     \item \code{Beta_raw_matrix}: k x V matrix of raw condition amplitudes
#'     \item \code{quality_metrics}: Diagnostics from the robust SVD
#'   }
#'   
#' @details This function implements Component 1, Step 4 of the M-HRF-LSS pipeline.
#'   It now delegates the SVD step to \code{extract_xi_beta_raw_svd_robust}, which
#'   adds condition number checks and fallback strategies. The first singular
#'   value and associated vectors are used to decompose gamma into Xi (HRF shape)
#'   and Beta (amplitude) components. Near-zero singular values are handled by
#'   setting the corresponding Xi and Beta values to zero. Quality metrics from
#'   the robust SVD can optionally be logged via \code{log_fun}.
#'   
#' @examples
#' \dontrun{
#' # Create example gamma coefficients
#' m <- 5   # manifold dimensions
#' k <- 3   # conditions
#' V <- 100 # voxels
#' 
#' # Gamma from GLM solve
#' gamma <- matrix(rnorm((k * m) * V), k * m, V)
#'
#' # Extract Xi and Beta
#' result <- extract_xi_beta_raw_svd_core(gamma, m, k)
#' # result$Xi_raw_matrix is m x V
#' # result$Beta_raw_matrix is k x V
#' }
#' 
#' @export
extract_xi_beta_raw_svd_core <- function(Gamma_coeffs_matrix,
                                        m_manifold_dim,
                                        k_conditions,
                                        n_jobs = 1,
                                        log_fun = NULL) {
  
  # Input validation
  if (!is.matrix(Gamma_coeffs_matrix)) {
    stop("Gamma_coeffs_matrix must be a matrix")
  }
  
  if (!is.numeric(m_manifold_dim) || length(m_manifold_dim) != 1 || 
      m_manifold_dim < 1 || m_manifold_dim != round(m_manifold_dim)) {
    stop("m_manifold_dim must be a positive integer")
  }
  
  if (!is.numeric(k_conditions) || length(k_conditions) != 1 || 
      k_conditions < 1 || k_conditions != round(k_conditions)) {
    stop("k_conditions must be a positive integer")
  }
  
  # Check dimensions
  expected_rows <- k_conditions * m_manifold_dim
  if (nrow(Gamma_coeffs_matrix) != expected_rows) {
    stop(sprintf(
      "Gamma_coeffs_matrix has %d rows but expected %d (k * m = %d * %d)",
      nrow(Gamma_coeffs_matrix), expected_rows, k_conditions, m_manifold_dim
    ))
  }
  
  V <- ncol(Gamma_coeffs_matrix)

  # Use robust SVD decomposition
  svd_result <- extract_xi_beta_raw_svd_robust(
    Gamma_coeffs_matrix = Gamma_coeffs_matrix,
    m_manifold_dim = m_manifold_dim,
    k_conditions = k_conditions,
    verbose_warnings = FALSE
  )

  Xi_raw_matrix <- svd_result$Xi_raw_matrix
  Beta_raw_matrix <- svd_result$Beta_raw_matrix

  # Optional logging of quality metrics
  if (!is.null(log_fun) && is.function(log_fun)) {
    gap_mean <- mean(svd_result$quality_metrics$singular_value_gaps, na.rm = TRUE)
    log_fun(sprintf("Mean singular value gap: %.4f", gap_mean))
  }

  # Return results with metrics
  list(
    Xi_raw_matrix = Xi_raw_matrix,
    Beta_raw_matrix = Beta_raw_matrix,
    quality_metrics = svd_result$quality_metrics
  )
}

#' Apply Intrinsic Identifiability (Core)
#'
#' Applies sign and scale constraints to ensure identifiability of HRF estimates.
#'
#' @param Xi_raw_matrix The m x V raw manifold coordinates from extract_xi_beta_raw_svd_core
#' @param Beta_raw_matrix The k x V raw condition amplitudes from extract_xi_beta_raw_svd_core
#' @param B_reconstructor_matrix The p x m HRF reconstructor matrix
#' @param h_ref_shape_vector The p x 1 reference HRF shape (e.g., canonical HRF)
#' @param ident_scale_method Scaling method: "l2_norm" (unit length), "max_abs_val" 
#'   (peak = 1), or "none"
#' @param ident_sign_method Sign alignment method: "canonical_correlation" (align with
#'   reference) or "data_fit_correlation" (requires additional data, not implemented)
#' @param zero_tol Numeric tolerance for treating a reconstructed HRF as zero
#'   when applying scale normalization. If the L2 norm (for "l2_norm") or
#'   maximum absolute value (for "max_abs_val") of a voxel's HRF is below this
#'   threshold, both \code{Xi_ident_matrix} and \code{Beta_ident_matrix} are set
#'   to zero for that voxel.
#' @param correlation_threshold Numeric threshold for canonical correlation
#'   (default 1e-3). When the absolute correlation between a voxel's HRF and
#'   the reference HRF falls below this threshold, the sign is determined
#'   using the RMS rule (sign of sum of HRF values) instead.
#' @param consistency_check Logical; if TRUE, the reconstructed HRF is
#'   re-projected after sign/scale alignment and verified against the
#'   canonical shape. Voxels failing the check are flipped to ensure
#'   consistent orientation.
#' @param n_jobs Number of parallel jobs for voxel processing (default 1).
#'   Only used when ident_sign_method = "data_fit_correlation".
#' @param verbose Logical; if TRUE, print progress messages for each voxel.
#'   Default FALSE to avoid message spam.
#'   
#' @return A list containing:
#'   \itemize{
#'     \item \code{Xi_ident_matrix}: m x V matrix of identifiability-constrained manifold coordinates
#'     \item \code{Beta_ident_matrix}: k x V matrix of identifiability-constrained condition amplitudes
#'   }
#'   
#' @details This function implements Component 1, Step 5 of the M-HRF-LSS pipeline.
#'   It ensures that HRF estimates are identifiable by:
#'   1. Aligning sign with a reference HRF (avoiding arbitrary sign flips)
#'   2. Normalizing scale consistently across voxels
#'   The chosen sign alignment method is reported via \code{message()}.
#'   If the absolute correlation with the canonical HRF is below
#'   \code{1e-3}, a warning is issued and the sign is determined using
#'   the RMS rule (sign of the sum of the HRF). When
#'   \code{consistency_check = TRUE}, each voxel's HRF is reprojected and the
#'   alignment with the canonical shape is reverified.
#'   The beta amplitudes are adjusted inversely to preserve the overall signal.
#'   Voxels whose reconstructed HRFs fall below \code{zero_tol} are zeroed in
#'   both returned matrices.
#'   
#' @examples
#' \dontrun{
#' # After SVD extraction
#' m <- 5
#' k <- 3
#' V <- 100
#' p <- 30
#' 
#' # Get raw Xi and Beta from SVD
#' svd_result <- extract_xi_beta_raw_svd_core(gamma, m, k)
#' 
#' # Apply identifiability with canonical HRF reference
#' h_canonical <- c(0, 0.8, 1, 0.7, 0.3, rep(0, p-5))  # simplified canonical
#' ident_result <- apply_intrinsic_identifiability_core(
#'   svd_result$Xi_raw_matrix,
#'   svd_result$Beta_raw_matrix,
#'   B_reconstructor,
#'   h_canonical
#' )
#' }
#' 
#' @export
apply_intrinsic_identifiability_core <- function(Xi_raw_matrix,
                                                Beta_raw_matrix,
                                                B_reconstructor_matrix,
                                                h_ref_shape_vector,
                                                ident_scale_method = "l2_norm",
                                                ident_sign_method = "canonical_correlation",
                                                zero_tol = 1e-8,
                                                correlation_threshold = 1e-3,
                                                Y_proj_matrix = NULL,
                                                X_condition_list_proj_matrices = NULL,
                                                consistency_check = FALSE,
                                                n_jobs = 1,
                                                verbose = FALSE) {

  # Input validation
  if (!is.matrix(Xi_raw_matrix)) {
    stop("Xi_raw_matrix must be a matrix")
  }
  
  if (!is.matrix(Beta_raw_matrix)) {
    stop("Beta_raw_matrix must be a matrix")
  }
  
  if (!is.matrix(B_reconstructor_matrix)) {
    stop("B_reconstructor_matrix must be a matrix")
  }
  
  if (!is.numeric(h_ref_shape_vector) || !is.vector(h_ref_shape_vector)) {
    stop("h_ref_shape_vector must be a numeric vector")
  }
  
  # Check dimensions
  m <- nrow(Xi_raw_matrix)
  V <- ncol(Xi_raw_matrix)
  k <- nrow(Beta_raw_matrix)
  p <- nrow(B_reconstructor_matrix)
  
  if (ncol(Beta_raw_matrix) != V) {
    stop("Xi_raw_matrix and Beta_raw_matrix must have the same number of columns")
  }
  
  if (ncol(B_reconstructor_matrix) != m) {
    stop("B_reconstructor_matrix must have m columns to match Xi dimension")
  }
  
  if (length(h_ref_shape_vector) != p) {
    stop("h_ref_shape_vector must have length p to match B_reconstructor rows")
  }
  
  # Validate method choices
  valid_scale_methods <- c("l2_norm", "max_abs_val", "none")
  if (!ident_scale_method %in% valid_scale_methods) {
    stop("ident_scale_method must be one of: ", paste(valid_scale_methods, collapse = ", "))
  }
  
  valid_sign_methods <- c("canonical_correlation", "data_fit_correlation")
  if (!ident_sign_method %in% valid_sign_methods) {
    stop("ident_sign_method must be one of: ", paste(valid_sign_methods, collapse = ", "))
  }

  message(sprintf("Using '%s' for sign alignment", ident_sign_method))
  
  # Use vectorized path for canonical correlation method
  if (ident_sign_method == "canonical_correlation") {
    # Efficient vectorized implementation
    result <- .apply_identifiability_vectorized(
      Xi_raw_matrix = Xi_raw_matrix,
      Beta_raw_matrix = Beta_raw_matrix,
      B_reconstructor_matrix = B_reconstructor_matrix,
      h_ref_shape_vector = h_ref_shape_vector,
      scale_method = ident_scale_method,
      correlation_threshold = correlation_threshold,
      zero_tol = zero_tol,
      consistency_check = consistency_check
    )
    
    return(result)
    
  } else if (ident_sign_method == "data_fit_correlation") {
    # Data fit method requires voxel-by-voxel processing
    if (is.null(Y_proj_matrix) || is.null(X_condition_list_proj_matrices)) {
      stop("data_fit_correlation method requires Y_proj_matrix and X_condition_list_proj_matrices")
    }
    
    # Compute reference manifold coordinates (unused but kept for compatibility)
    xi_ref_coord <- MASS::ginv(B_reconstructor_matrix) %*% h_ref_shape_vector
    
    # Set up configuration for voxel processing
    config <- list(
      sign_method = ident_sign_method,
      scale_method = ident_scale_method,
      zero_tol = zero_tol,
      correlation_threshold = correlation_threshold,
      consistency_check = consistency_check
    )
    
    # Process each voxel using the helper function
    voxel_fun <- function(vx) {
      .process_voxel_identifiability(
        vx = vx,
        Xi_raw = Xi_raw_matrix,
        Beta_raw = Beta_raw_matrix,
        B_reconstructor = B_reconstructor_matrix,
        h_ref = h_ref_shape_vector,
        config = config,
        Y_proj = Y_proj_matrix,
        X_list = X_condition_list_proj_matrices,
        verbose = verbose
      )
    }
    
    # Process voxels in parallel if requested
    res_list <- .parallel_lapply(seq_len(V), voxel_fun, n_jobs)
    
    # Assemble results
    Xi_ident_matrix <- do.call(cbind, lapply(res_list, "[[", "xi"))
    Beta_ident_matrix <- do.call(cbind, lapply(res_list, "[[", "beta"))
    
    # Return results
    list(
      Xi_ident_matrix = Xi_ident_matrix,
      Beta_ident_matrix = Beta_ident_matrix
    )
  }
}
