# Core LSS (Least Squares Separate) Functions (Component 3)
# Implementation of MHRF-CORE-LSS-01 through MHRF-CORE-LSS-04

#' Prepare LSS Fixed Components (Core)
#'
#' Precomputes fixed regression components for efficient Least Squares Separate 
#' (LSS) estimation using the Woodbury matrix identity.
#'
#' @param A_lss_fixed_matrix An n x q_lss matrix of fixed regressors (e.g., 
#'   intercept, drift terms, nuisance regressors), where n is timepoints and 
#'   q_lss is number of fixed regressors
#' @param intercept_col_index_in_Alss Integer index of the intercept column in 
#'   A_lss_fixed_matrix, or NULL if no intercept
#' @param lambda_ridge_Alss Ridge penalty parameter for regularizing the fixed 
#'   regressor inversion (scalar, typically small like 1e-6)
#'   
#' @return A list containing:
#'   \itemize{
#'     \item \code{P_lss_matrix}: q_lss x n matrix equal to 
#'       (A'A + lambda*I)^(-1) * A'
#'     \item \code{p_lss_vector}: n x 1 vector for intercept projection. If 
#'       intercept_col_index_in_Alss is provided, this is the corresponding 
#'       row of the Moore-Penrose pseudoinverse of A. Otherwise, a zero vector.
#'   }
#'   
#' @details This function implements Component 3, Step 1 of the M-HRF-LSS pipeline.
#'   It precomputes matrices needed for the Woodbury identity-based efficient
#'   computation of single-trial betas. The P_lss matrix projects data onto the
#'   fixed regressor space, while p_lss_vector handles the special case of the
#'   intercept term in the LSS formulation.
#'   
#' @examples
#' \dontrun{
#' # Create example fixed regressors
#' n <- 200  # timepoints
#' 
#' # Intercept + linear drift + motion parameters
#' A_fixed <- cbind(
#'   1,                    # intercept
#'   seq_len(n) / n,       # linear drift
#'   rnorm(n),             # motion param 1
#'   rnorm(n)              # motion param 2
#' )
#' 
#' # Prepare LSS components
#' lss_components <- prepare_lss_fixed_components_core(
#'   A_fixed, 
#'   intercept_col_index_in_Alss = 1,
#'   lambda_ridge_Alss = 1e-6
#' )
#' }
#' 
#' @export
prepare_lss_fixed_components_core <- function(A_lss_fixed_matrix,
                                            intercept_col_index_in_Alss = NULL,
                                            lambda_ridge_Alss = 0) {
  
  # Input validation
  if (!is.matrix(A_lss_fixed_matrix)) {
    stop("A_lss_fixed_matrix must be a matrix")
  }
  
  n <- nrow(A_lss_fixed_matrix)
  q_lss <- ncol(A_lss_fixed_matrix)
  
  if (n < 2) {
    stop("A_lss_fixed_matrix must have at least 2 rows (timepoints)")
  }
  
  if (q_lss < 1) {
    stop("A_lss_fixed_matrix must have at least 1 column")
  }
  
  if (q_lss >= n) {
    stop(sprintf(
      "A_lss_fixed_matrix has too many columns (%d) relative to rows (%d). Must have q_lss < n.",
      q_lss, n
    ))
  }
  
  if (!is.null(intercept_col_index_in_Alss)) {
    if (!is.numeric(intercept_col_index_in_Alss) || 
        length(intercept_col_index_in_Alss) != 1 ||
        intercept_col_index_in_Alss < 1 || 
        intercept_col_index_in_Alss > q_lss ||
        intercept_col_index_in_Alss != round(intercept_col_index_in_Alss)) {
      stop(sprintf(
        "intercept_col_index_in_Alss must be an integer between 1 and %d", 
        q_lss
      ))
    }
  }
  
  lambda_ridge_Alss <- .validate_and_standardize_lambda(lambda_ridge_Alss,
                                                        "lambda_ridge_Alss")
  
  # Step 1: Compute A'A
  AtA <- crossprod(A_lss_fixed_matrix)  # q_lss x q_lss
  
  # Step 2: Add ridge penalty to diagonal
  AtA_reg <- AtA + lambda_ridge_Alss * diag(q_lss)
  
  # Step 3: Check condition number and add jitter if needed
  # This helps with near-singular matrices
  eigvals <- eigen(AtA_reg, symmetric = TRUE, only.values = TRUE)$values
  min_eigval <- min(eigvals)
  max_eigval <- max(eigvals)
  condition_number <- max_eigval / max(min_eigval, .Machine$double.eps)
  
  if (condition_number > 1e8) {
    # Add small jitter to improve conditioning
    jitter <- 1e-6 * median(diag(AtA_reg))
    AtA_reg <- AtA_reg + jitter * diag(q_lss)
    warning(sprintf(
      "Fixed regressor matrix is poorly conditioned (condition number = %.2e). Added jitter = %.2e to diagonal.",
      condition_number, jitter
    ))
  }
  
  # Step 4: Compute P_lss = (A'A + lambda*I)^(-1) * A'
  # Use Cholesky decomposition for stable inversion
  AtA_reg_inv <- chol2inv(chol(AtA_reg))

  # Then multiply by A'
  P_lss_matrix <- AtA_reg_inv %*% t(A_lss_fixed_matrix)  # q_lss x n
  
  # Step 5: Compute p_lss_vector for intercept handling
  if (!is.null(intercept_col_index_in_Alss)) {
    # Compute pseudoinverse of A_lss_fixed_matrix
    # Using the already computed components for efficiency
    # A+ = (A'A)^(-1) * A' when A has full column rank
    # We already have this as P_lss_matrix

    # Extract the row corresponding to the intercept
    p_lss_vector <- P_lss_matrix[intercept_col_index_in_Alss, , drop = TRUE]
  } else {
    # No intercept specified, use zero vector
    p_lss_vector <- rep(0, n)
  }
  
  # Return precomputed components
  list(
    P_lss_matrix = P_lss_matrix,
    p_lss_vector = p_lss_vector
  )
}

#' Reconstruct HRF Shapes from Manifold Coordinates (Core)
#'
#' Reconstructs voxel-specific HRF shapes by multiplying the HRF basis 
#' reconstructor matrix with spatially smoothed manifold coordinates.
#'
#' @param B_reconstructor_matrix The p x m HRF reconstructor matrix from 
#'   Component 0, where p is HRF length and m is manifold dimensionality
#' @param Xi_smoothed_matrix The m x V matrix of spatially smoothed manifold 
#'   coordinates from Component 2, where V is number of voxels
#'   
#' @return H_shapes_allvox_matrix A p x V matrix where each column is the 
#'   reconstructed HRF shape for a voxel
#'   
#' @details This function implements Component 3, Step 2 of the M-HRF-LSS pipeline.
#'   It transforms the low-dimensional manifold representation back to the full
#'   HRF shape space. The resulting HRF shapes are used to create voxel-specific
#'   regressors for the LSS estimation.
#'   
#' @examples
#' \dontrun{
#' # Example with synthetic data
#' p <- 30   # HRF length
#' m <- 5    # manifold dimensions
#' V <- 100  # voxels
#' 
#' # Example reconstructor (would come from Component 0)
#' B_reconstructor <- matrix(rnorm(p * m), p, m)
#' 
#' # Example smoothed manifold coordinates (would come from Component 2)
#' Xi_smoothed <- matrix(rnorm(m * V), m, V)
#' 
#' # Reconstruct HRF shapes
#' H_shapes <- reconstruct_hrf_shapes_core(B_reconstructor, Xi_smoothed)
#' # Result is 30 x 100 matrix (HRF samples x voxels)
#' }
#' 
#' @export
reconstruct_hrf_shapes_core <- function(B_reconstructor_matrix,
                                      Xi_smoothed_matrix) {
  
  # Input validation
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
  
  if (p < 2) {
    stop("B_reconstructor_matrix must have at least 2 rows (HRF time points)")
  }
  
  if (m_B < 1) {
    stop("B_reconstructor_matrix must have at least 1 column (manifold dimension)")
  }
  
  if (V < 1) {
    stop("Xi_smoothed_matrix must have at least 1 column (voxel)")
  }
  
  # Reconstruct HRF shapes
  # H = B * Xi
  # where B is p x m and Xi is m x V
  # resulting in H that is p x V
  H_shapes_allvox_matrix <- B_reconstructor_matrix %*% Xi_smoothed_matrix
  
  # Ensure the result is a regular matrix (not Matrix class)
  if (inherits(H_shapes_allvox_matrix, "Matrix")) {
    H_shapes_allvox_matrix <- as.matrix(H_shapes_allvox_matrix)
  }
  
  return(H_shapes_allvox_matrix)
}

#' Compute trial-wise betas using the Woodbury identity
#'
#' Internal helper used by LSS functions to compute single-trial betas
#' given the voxel-specific trial regressor matrix.
#'
#' @param C_v Matrix of convolved trial regressors (n x T)
#' @param Y_v Numeric vector of length n with projected voxel data
#' @param A_fixed Fixed regressors matrix used during projection (n x q)
#' @param P_lss Precomputed matrix from `prepare_lss_fixed_components_core`
#'              with dimensions q x n
#' @param p_lss Precomputed intercept projection vector of length n
#'
#' @return Numeric vector of length T containing LSS beta estimates
#' @keywords internal
.compute_lss_betas <- function(C_v, Y_v, A_fixed, P_lss, p_lss) {
  U_v <- P_lss %*% C_v
  V_regressors_v <- C_v - A_fixed %*% U_v
  pc_v_row <- as.vector(crossprod(p_lss, C_v))
  cv_v_row <- colSums(V_regressors_v * V_regressors_v)
  alpha_v_row <- (1 - pc_v_row) / pmax(cv_v_row, .Machine$double.eps)
  S_effective_regressors_v <- sweep(V_regressors_v, MARGIN = 2,
                                    STATS = alpha_v_row, FUN = "*")
  S_effective_regressors_v <- sweep(S_effective_regressors_v, MARGIN = 1,
                                    STATS = p_lss, FUN = "+")
  as.vector(crossprod(S_effective_regressors_v, Y_v))
}

#' Run LSS for Single Voxel (Core)

#' Run LSS for Single Voxel (Core)
#'
#' Performs Least Squares Separate (LSS) estimation for a single voxel using
#' the Woodbury matrix identity for computational efficiency.
#'
#' @param Y_proj_voxel_vector An n x 1 vector of projected BOLD data for one voxel
#' @param X_trial_onset_list_of_matrices A list of T matrices, where each is an 
#'   n x p Toeplitz design matrix for a single trial
#' @param H_shape_voxel_vector A p x 1 vector representing the HRF shape for 
#'   this voxel
#' @param A_lss_fixed_matrix The n x q_lss matrix of fixed regressors
#' @param P_lss_matrix The q_lss x n precomputed projection matrix from 
#'   prepare_lss_fixed_components_core
#' @param p_lss_vector The n x 1 precomputed intercept projection vector from 
#'   prepare_lss_fixed_components_core
#'   
#' @return beta_trial_voxel_vector A T x 1 vector of trial-wise beta estimates
#'   
#' @details This function implements Component 3, Step 3 of the M-HRF-LSS pipeline.
#'   For each trial, it uses the Woodbury matrix identity to efficiently compute
#'   the beta estimate when that trial's regressor is added to the fixed model.
#'   This avoids repeated matrix inversions and makes single-trial estimation
#'   computationally feasible.
#'   
#'   The Woodbury identity allows us to update the inverse when adding a single
#'   column c to the design matrix:
#'   (A'A + cc')^(-1) = (A'A)^(-1) - (A'A)^(-1)cc'(A'A)^(-1) / (1 + c'(A'A)^(-1)c)
#'   
#' @examples
#' \dontrun{
#' # Setup for single voxel
#' n <- 200  # timepoints
#' p <- 30   # HRF length
#' T <- 50   # trials
#' 
#' # Data for one voxel
#' Y_voxel <- rnorm(n)
#' 
#' # Trial onsets (simplified)
#' X_trials <- lapply(1:T, function(t) {
#'   X <- matrix(0, n, p)
#'   # Put onset at different times
#'   onset_time <- 10 + (t-1) * 3
#'   if (onset_time + p <= n) {
#'     X[onset_time:(onset_time+p-1), ] <- diag(p)
#'   }
#'   X
#' })
#' 
#' # HRF shape for this voxel
#' H_voxel <- dgamma(0:(p-1), shape = 6, rate = 1)
#' H_voxel <- H_voxel / max(H_voxel)
#' 
#' # Fixed regressors
#' A_fixed <- cbind(1, seq_len(n)/n)
#' lss_prep <- prepare_lss_fixed_components_core(A_fixed, 1, 1e-6)
#' 
#' # Run LSS
#' betas <- run_lss_for_voxel_core(
#'   Y_voxel, X_trials, H_voxel, A_fixed,
#'   lss_prep$P_lss_matrix, lss_prep$p_lss_vector
#' )
#' }
#' 
#' @export
run_lss_for_voxel_core <- function(Y_proj_voxel_vector,
                                  X_trial_onset_list_of_matrices,
                                  H_shape_voxel_vector,
                                  A_lss_fixed_matrix,
                                  P_lss_matrix,
                                  p_lss_vector) {
  
  # Accept single-column matrices for vector input
  if (is.matrix(Y_proj_voxel_vector) && ncol(Y_proj_voxel_vector) == 1) {
    Y_proj_voxel_vector <- drop(Y_proj_voxel_vector)
  }
  # Input validation
  if (!is.numeric(Y_proj_voxel_vector) || !is.vector(Y_proj_voxel_vector)) {
    stop("Y_proj_voxel_vector must be a numeric vector")
  }
  
  n <- length(Y_proj_voxel_vector)
  
  if (!is.list(X_trial_onset_list_of_matrices)) {
    stop("X_trial_onset_list_of_matrices must be a list")
  }
  
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  if (T_trials < 1) {
    stop("X_trial_onset_list_of_matrices must contain at least one trial")
  }
  
  if (!is.numeric(H_shape_voxel_vector) || !is.vector(H_shape_voxel_vector)) {
    stop("H_shape_voxel_vector must be a numeric vector")
  }
  
  p <- length(H_shape_voxel_vector)
  
  if (!is.matrix(A_lss_fixed_matrix)) {
    stop("A_lss_fixed_matrix must be a matrix")
  }
  
  if (nrow(A_lss_fixed_matrix) != n) {
    stop("A_lss_fixed_matrix must have n rows to match Y_proj_voxel_vector length")
  }
  
  q_lss <- ncol(A_lss_fixed_matrix)
  
  if (!is.matrix(P_lss_matrix)) {
    stop("P_lss_matrix must be a matrix")
  }
  
  if (nrow(P_lss_matrix) != q_lss || ncol(P_lss_matrix) != n) {
    stop(sprintf("P_lss_matrix must be %d x %d", q_lss, n))
  }
  
  if (!is.numeric(p_lss_vector) || !is.vector(p_lss_vector) || length(p_lss_vector) != n) {
    stop("p_lss_vector must be a numeric vector of length n")
  }
  
  # Validate each trial matrix
  for (t in 1:T_trials) {
    if (!is.matrix(X_trial_onset_list_of_matrices[[t]])) {
      stop(sprintf("X_trial_onset_list_of_matrices[[%d]] must be a matrix", t))
    }
    if (nrow(X_trial_onset_list_of_matrices[[t]]) != n) {
      stop(sprintf("X_trial_onset_list_of_matrices[[%d]] must have %d rows", t, n))
    }
    if (ncol(X_trial_onset_list_of_matrices[[t]]) != p) {
      stop(sprintf("X_trial_onset_list_of_matrices[[%d]] must have %d columns", t, p))
    }
  }
  
  # Step 1: Construct C_v matrix (n x T)
  # Each column is the convolution of trial onset with voxel-specific HRF
  C_v <- matrix(0, nrow = n, ncol = T_trials)
  
  for (t in 1:T_trials) {
    # Multiply trial design matrix by HRF shape
    # X_t is n x p, H is p x 1, result is n x 1
    C_v[, t] <- X_trial_onset_list_of_matrices[[t]] %*% H_shape_voxel_vector
  }
  
  # Step 2: Woodbury LSS computation via helper
  beta_trial_voxel_vector <- .compute_lss_betas(
    C_v,
    Y_proj_voxel_vector,
    A_lss_fixed_matrix,
    P_lss_matrix,
    p_lss_vector
  )

  return(beta_trial_voxel_vector)
}

#'
#' @param Y_proj_matrix n x V projected BOLD data matrix, where n is
#'   the number of timepoints and V is the number of voxels.
#' @param X_trial_onset_list_of_matrices List of length T, each an n x p matrix
#'   of trial-specific regressors, where p is the HRF length.
#' @param H_shapes_allvox_matrix p x V matrix of HRF shapes for all voxels.
#' @param A_lss_fixed_matrix n x q matrix of fixed regressors to residualize.
#' @param P_lss_matrix n x n residualizing projection matrix from fixed regressors.
#' @param p_lss_vector Numeric vector of length n for residualization of fixed effects.
#' @param ram_heuristic_GB_for_Rt Numeric scalar indicating RAM limit in gigabytes
#'   for precomputing components (default 1.0).
#' @param use_fmrireg Logical; if \code{TRUE}, uses \pkg{fmrireg} internals
#'   for model fitting (default \code{TRUE}).
#' @param n_jobs Number of parallel jobs for voxel processing (default 1).
#' @return A T x V matrix of trial-level beta estimates, where T is the number
#'   of trials and V is the number of voxels.
#' @export
run_lss_voxel_loop_core <- function(Y_proj_matrix,
                                   X_trial_onset_list_of_matrices,
                                   H_shapes_allvox_matrix,
                                   A_lss_fixed_matrix,
                                   P_lss_matrix,
                                   p_lss_vector,
                                   ram_heuristic_GB_for_Rt = 1.0,
                                   use_fmrireg = TRUE,
                                   n_jobs = 1) {
  
  # Input validation
  if (!is.matrix(Y_proj_matrix)) {
    stop("Y_proj_matrix must be a matrix")
  }
  
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  
  if (!is.list(X_trial_onset_list_of_matrices)) {
    stop("X_trial_onset_list_of_matrices must be a list")
  }
  
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  if (T_trials < 1) {
    stop("X_trial_onset_list_of_matrices must contain at least one trial")
  }
  
  if (!is.matrix(H_shapes_allvox_matrix)) {
    stop("H_shapes_allvox_matrix must be a matrix")
  }
  
  p <- nrow(H_shapes_allvox_matrix)
  
  if (ncol(H_shapes_allvox_matrix) != V) {
    stop("H_shapes_allvox_matrix must have V columns to match Y_proj_matrix")
  }
  
  if (!is.matrix(A_lss_fixed_matrix)) {
    stop("A_lss_fixed_matrix must be a matrix")
  }
  
  if (nrow(A_lss_fixed_matrix) != n) {
    stop("A_lss_fixed_matrix must have n rows to match Y_proj_matrix")
  }
  
  if (!is.matrix(P_lss_matrix)) {
    stop("P_lss_matrix must be a matrix")
  }
  
  if (!is.numeric(p_lss_vector) || !is.vector(p_lss_vector) || length(p_lss_vector) != n) {
    stop("p_lss_vector must be a numeric vector of length n")
  }
  
  if (!is.numeric(ram_heuristic_GB_for_Rt) || length(ram_heuristic_GB_for_Rt) != 1 || 
      ram_heuristic_GB_for_Rt < 0) {
    stop("ram_heuristic_GB_for_Rt must be a non-negative scalar")
  }
  
  # Validate trial matrices
  for (t in 1:T_trials) {
    if (!is.matrix(X_trial_onset_list_of_matrices[[t]])) {
      stop(sprintf("X_trial_onset_list_of_matrices[[%d]] must be a matrix", t))
    }
    if (nrow(X_trial_onset_list_of_matrices[[t]]) != n ||
        ncol(X_trial_onset_list_of_matrices[[t]]) != p) {
      stop(sprintf("X_trial_onset_list_of_matrices[[%d]] must be %d x %d", t, n, p))
    }
  }
  
  # Initialize output matrix
  Beta_trial_allvox_matrix <- matrix(0, nrow = T_trials, ncol = V)
  
  # Try to use fmrireg implementation if available and requested
  if (use_fmrireg && requireNamespace("fmrireg", quietly = TRUE)) {
    # Use fmrireg's lss_compute_r function which works on pre-projected data
    # This is the same algorithm but potentially more optimized
    
    # We need to compute Q_dmat_ran which is the projected trial regressors
    # Since Y_proj_matrix is already projected by P_confound = I - A(A'A)^{-1}A'
    # We need to apply the same projection to the trial regressors
    
    # First compute the projection matrix if not provided
    # P_confound = I - A_lss_fixed %*% ginv(A_lss_fixed)
    P_confound <- diag(n) - A_lss_fixed_matrix %*% MASS::ginv(A_lss_fixed_matrix)
    
    voxel_fun <- function(v) {
      dmat_ran <- matrix(0, n, T_trials)
      for (t in 1:T_trials) {
        dmat_ran[, t] <- X_trial_onset_list_of_matrices[[t]] %*% H_shapes_allvox_matrix[, v]
      }
      Q_dmat_ran <- P_confound %*% dmat_ran
      residual_data <- Y_proj_matrix[, v, drop = FALSE]
      as.vector(fmrireg:::lss_compute_r(Q_dmat_ran, residual_data))
    }

    res_list <- .parallel_lapply(seq_len(V), voxel_fun, n_jobs)
    Beta_trial_allvox_matrix <- do.call(cbind, res_list)

    return(Beta_trial_allvox_matrix)
  }
  
  # Memory heuristic: Check if we can precompute all R_t matrices
  # Each R_t is n x V, we have T of them
  # Memory in bytes: T * n * V * 8 (assuming double precision)
  estimated_memory_GB <- (T_trials * n * V * 8) / (1024^3)
  
  precompute_R_t <- (estimated_memory_GB < ram_heuristic_GB_for_Rt)
  
  if (precompute_R_t) {
    # Precompute all R_t matrices for efficiency
    # R_t = X_t %*% H_shapes for all trials
    R_t_allvox_list <- vector("list", T_trials)
    
    for (t in 1:T_trials) {
      # X_t is n x p, H_shapes is p x V, result is n x V
      R_t_allvox_list[[t]] <- X_trial_onset_list_of_matrices[[t]] %*% H_shapes_allvox_matrix
    }
  }
  
  voxel_fun <- function(v) {
    Y_voxel <- Y_proj_matrix[, v]
    H_voxel <- H_shapes_allvox_matrix[, v]

    if (precompute_R_t) {
      C_v <- matrix(0, nrow = n, ncol = T_trials)
      for (t in 1:T_trials) {
        C_v[, t] <- R_t_allvox_list[[t]][, v]
      }

      .compute_lss_betas(C_v, Y_voxel, A_lss_fixed_matrix,
                         P_lss_matrix, p_lss_vector)
    } else {
      run_lss_for_voxel_core(
        Y_voxel,
        X_trial_onset_list_of_matrices,
        H_voxel,
        A_lss_fixed_matrix,
        P_lss_matrix,
        p_lss_vector
      )
    }
  }

  res_list <- .parallel_lapply(seq_len(V), voxel_fun, n_jobs)
  Beta_trial_allvox_matrix <- do.call(cbind, res_list)
  
  return(Beta_trial_allvox_matrix)
}
