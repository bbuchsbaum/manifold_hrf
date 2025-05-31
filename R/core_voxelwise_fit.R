# Core Voxel-wise Manifold Fit Functions (Component 1)
# Implementation of MHRF-CORE-VOXFIT-01 through MHRF-CORE-VOXFIT-05

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
#'   It uses QR decomposition to efficiently project out confound regressors from
#'   both the data and design matrices. If Z_confounds_matrix is NULL, the 
#'   original matrices are returned unchanged.
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
  
  if (!is.list(X_list_of_matrices)) {
    stop("X_list_of_matrices must be a list")
  }
  
  n <- nrow(Y_data_matrix)
  
  # Check that all X matrices have the same number of rows as Y
  for (i in seq_along(X_list_of_matrices)) {
    if (!is.matrix(X_list_of_matrices[[i]])) {
      stop(sprintf("X_list_of_matrices[[%d]] must be a matrix", i))
    }
    if (nrow(X_list_of_matrices[[i]]) != n) {
      stop(sprintf("X_list_of_matrices[[%d]] must have %d rows to match Y_data_matrix", i, n))
    }
  }
  
  # If no confounds, return original matrices
  if (is.null(Z_confounds_matrix)) {
    return(list(
      Y_proj_matrix = Y_data_matrix,
      X_list_proj_matrices = X_list_of_matrices
    ))
  }
  
  # Validate confounds matrix
  if (!is.matrix(Z_confounds_matrix)) {
    stop("Z_confounds_matrix must be a matrix or NULL")
  }
  
  if (nrow(Z_confounds_matrix) != n) {
    stop("Z_confounds_matrix must have the same number of rows as Y_data_matrix")
  }
  
  # Check for rank deficiency in confounds
  if (ncol(Z_confounds_matrix) >= n) {
    stop("Z_confounds_matrix has too many columns (must be less than number of timepoints)")
  }
  
  # Step 1: Compute QR decomposition of confounds matrix
  # Q_Z is an orthonormal basis for the column space of Z
  qr_Z <- qr(Z_confounds_matrix, LAPACK = TRUE)
  Q_Z <- qr.Q(qr_Z)
  
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
  
  if (!is.numeric(lambda_gamma) || length(lambda_gamma) != 1 || lambda_gamma < 0) {
    stop("lambda_gamma must be a non-negative scalar")
  }
  
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
  
  # Step 5: Compute X'Y
  XtY <- crossprod(X_tilde, Y_proj_matrix)  # (k*m) x V
  
  # Step 6: Solve for gamma coefficients
  # Gamma = (X'X + lambda*I)^(-1) * X'Y
  Gamma_coeffs_matrix <- solve(XtX_reg, XtY)
  
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
#' 
#' @return A list containing:
#'   \itemize{
#'     \item \code{Xi_raw_matrix}: m x V matrix of raw manifold coordinates
#'     \item \code{Beta_raw_matrix}: k x V matrix of raw condition amplitudes
#'   }
#'   
#' @details This function implements Component 1, Step 4 of the M-HRF-LSS pipeline.
#'   For each voxel, it reshapes the gamma coefficients into an m x k matrix and
#'   performs SVD to extract the dominant pattern. The first singular value and
#'   vectors are used to decompose gamma into Xi (HRF shape) and Beta (amplitude)
#'   components. Near-zero singular values are handled by setting the corresponding
#'   Xi and Beta values to zero.
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
                                        k_conditions) {
  
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
  
  # Initialize output matrices
  Xi_raw_matrix <- matrix(0, nrow = m_manifold_dim, ncol = V)
  Beta_raw_matrix <- matrix(0, nrow = k_conditions, ncol = V)
  
  # Threshold for near-zero singular values
  svd_threshold <- sqrt(.Machine$double.eps)
  
  # Loop over voxels
  for (vx in 1:V) {
    # Extract gamma coefficients for this voxel
    gamma_vx <- Gamma_coeffs_matrix[, vx]
    
    # Reshape into m x k matrix
    # The gamma vector is organized as:
    # [cond1_dim1, cond1_dim2, ..., cond1_dimM, cond2_dim1, ..., cond2_dimM, ...]
    # We want to reshape so rows are manifold dimensions and columns are conditions
    G_vx <- matrix(gamma_vx, nrow = m_manifold_dim, ncol = k_conditions)
    
    # Perform SVD
    svd_result <- svd(G_vx)
    
    # Check if first singular value is above threshold
    if (svd_result$d[1] > svd_threshold) {
      # Extract first singular components
      # Split singular value between u and v using square root
      sqrt_s1 <- sqrt(svd_result$d[1])
      
      # Xi gets u * sqrt(s)
      Xi_raw_matrix[, vx] <- svd_result$u[, 1] * sqrt_s1
      
      # Beta gets v * sqrt(s)
      Beta_raw_matrix[, vx] <- svd_result$v[, 1] * sqrt_s1
    }
    # If singular value is near zero, Xi and Beta remain zero for this voxel
  }
  
  # Return results
  list(
    Xi_raw_matrix = Xi_raw_matrix,
    Beta_raw_matrix = Beta_raw_matrix
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
#'   The beta amplitudes are adjusted inversely to preserve the overall signal.
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
                                                ident_sign_method = "canonical_correlation") {
  
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
  
  if (ident_sign_method == "data_fit_correlation") {
    stop("data_fit_correlation method requires additional data inputs and is not yet implemented")
  }
  
  # Compute reference manifold coordinates
  # xi_ref = (B_reconstructor)^+ * h_ref
  # Using Moore-Penrose pseudoinverse
  xi_ref_coord <- MASS::ginv(B_reconstructor_matrix) %*% h_ref_shape_vector
  
  # Initialize output matrices
  Xi_ident_matrix <- matrix(0, nrow = m, ncol = V)
  Beta_ident_matrix <- matrix(0, nrow = k, ncol = V)
  
  # Process each voxel
  for (vx in 1:V) {
    xi_vx <- Xi_raw_matrix[, vx]
    beta_vx <- Beta_raw_matrix[, vx]
    
    # Skip if xi is zero (no signal)
    if (all(abs(xi_vx) < .Machine$double.eps)) {
      Xi_ident_matrix[, vx] <- xi_vx
      Beta_ident_matrix[, vx] <- beta_vx
      next
    }
    
    # Step 1: Sign alignment
    if (ident_sign_method == "canonical_correlation") {
      # Align sign based on correlation with reference
      sgn <- sign(sum(xi_vx * xi_ref_coord))
      # Handle case where correlation is exactly zero
      if (sgn == 0) sgn <- 1
    }
    
    # Apply sign
    xi_vx_signed <- xi_vx * sgn
    beta_vx_signed <- beta_vx * sgn
    
    # Reconstruct HRF for scaling
    reconstructed_hrf_vx <- B_reconstructor_matrix %*% xi_vx_signed
    
    # Step 2: Scale normalization
    if (ident_scale_method == "l2_norm") {
      # Scale by L2 norm
      l2_norm <- sqrt(sum(reconstructed_hrf_vx^2))
      scl <- 1 / max(l2_norm, .Machine$double.eps)
    } else if (ident_scale_method == "max_abs_val") {
      # Scale by maximum absolute value
      max_abs <- max(abs(reconstructed_hrf_vx))
      scl <- 1 / max(max_abs, .Machine$double.eps)
    } else {  # "none"
      scl <- 1
    }
    
    # Apply scaling
    # Xi is scaled up, Beta is scaled down to preserve signal
    Xi_ident_matrix[, vx] <- xi_vx_signed * scl
    Beta_ident_matrix[, vx] <- beta_vx_signed / scl
    
    # If scaling was based on machine epsilon (essentially zero HRF), set Beta to zero
    if (scl > 1 / sqrt(.Machine$double.eps)) {
      Beta_ident_matrix[, vx] <- 0
    }
  }
  
  # Return results
  list(
    Xi_ident_matrix = Xi_ident_matrix,
    Beta_ident_matrix = Beta_ident_matrix
  )
}