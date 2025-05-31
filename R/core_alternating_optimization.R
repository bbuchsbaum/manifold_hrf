# Core Alternating Optimization Functions (Component 4)
# Implementation of MHRF-CORE-ALTOPT-01

#' Estimate Final Condition Betas (Core)
#'
#' Re-estimates condition-level beta coefficients using the spatially smoothed
#' HRF shapes from the manifold estimation pipeline.
#'
#' @param Y_proj_matrix The n x V confound-projected data matrix, where n is 
#'   timepoints and V is number of voxels
#' @param X_condition_list_proj_matrices List of k projected design matrices 
#'   (each n x p), where k is number of conditions and p is HRF length
#' @param H_shapes_allvox_matrix The p x V matrix of smoothed voxel-specific 
#'   HRF shapes from Component 3
#' @param lambda_beta_final Ridge penalty parameter for final beta estimation 
#'   (scalar, typically small like 0.01)
#' @param control_alt_list List with control parameters:
#'   \itemize{
#'     \item \code{max_iter}: Maximum number of iterations (default 1 for MVP)
#'     \item \code{rel_change_tol}: Relative change tolerance for convergence 
#'       (default 1e-4)
#'   }
#'   
#' @return Beta_condition_final_matrix A k x V matrix of final condition-level 
#'   beta estimates
#'   
#' @details This function implements Component 4 of the M-HRF-LSS pipeline.
#'   It re-estimates condition-level amplitudes using the final smoothed HRF
#'   shapes. For each voxel, it constructs condition-specific regressors by
#'   convolving the design matrices with the voxel's HRF, then solves a
#'   ridge regression problem. The MVP version does a single pass (max_iter=1),
#'   but the framework supports iterative refinement between HRFs and betas.
#'   
#' @examples
#' \dontrun{
#' # Setup
#' n <- 200  # timepoints
#' p <- 30   # HRF length
#' k <- 3    # conditions
#' V <- 100  # voxels
#' 
#' # Projected data
#' Y_proj <- matrix(rnorm(n * V), n, V)
#' 
#' # Condition design matrices
#' X_cond_list <- lapply(1:k, function(c) {
#'   # Simple block design for each condition
#'   X <- matrix(0, n, p)
#'   # Add some events
#'   onsets <- seq(10 + (c-1)*20, n-p, by = 60)
#'   for (onset in onsets) {
#'     X[onset:(onset+p-1), ] <- diag(p)
#'   }
#'   X
#' })
#' 
#' # HRF shapes (from previous components)
#' H_shapes <- matrix(rnorm(p * V), p, V)
#' 
#' # Estimate final betas
#' Beta_final <- estimate_final_condition_betas_core(
#'   Y_proj, X_cond_list, H_shapes,
#'   lambda_beta_final = 0.01,
#'   control_alt_list = list(max_iter = 1)
#' )
#' }
#' 
#' @export
estimate_final_condition_betas_core <- function(Y_proj_matrix,
                                              X_condition_list_proj_matrices,
                                              H_shapes_allvox_matrix,
                                              lambda_beta_final = 0.01,
                                              control_alt_list = list(max_iter = 1, 
                                                                    rel_change_tol = 1e-4)) {
  
  # Input validation
  if (!is.matrix(Y_proj_matrix)) {
    stop("Y_proj_matrix must be a matrix")
  }
  
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  
  if (!is.list(X_condition_list_proj_matrices)) {
    stop("X_condition_list_proj_matrices must be a list")
  }
  
  k <- length(X_condition_list_proj_matrices)
  
  if (k < 1) {
    stop("X_condition_list_proj_matrices must contain at least one condition")
  }
  
  if (!is.matrix(H_shapes_allvox_matrix)) {
    stop("H_shapes_allvox_matrix must be a matrix")
  }
  
  p <- nrow(H_shapes_allvox_matrix)
  
  if (ncol(H_shapes_allvox_matrix) != V) {
    stop("H_shapes_allvox_matrix must have V columns to match Y_proj_matrix")
  }
  
  if (!is.numeric(lambda_beta_final) || length(lambda_beta_final) != 1 || 
      lambda_beta_final < 0) {
    stop("lambda_beta_final must be a non-negative scalar")
  }
  
  # Validate control parameters
  if (!is.list(control_alt_list)) {
    stop("control_alt_list must be a list")
  }
  
  # Set defaults for control parameters
  if (is.null(control_alt_list$max_iter)) {
    control_alt_list$max_iter <- 1
  }
  if (is.null(control_alt_list$rel_change_tol)) {
    control_alt_list$rel_change_tol <- 1e-4
  }
  
  max_iter <- control_alt_list$max_iter
  rel_change_tol <- control_alt_list$rel_change_tol
  
  if (!is.numeric(max_iter) || length(max_iter) != 1 || max_iter < 1 || 
      max_iter != round(max_iter)) {
    stop("max_iter must be a positive integer")
  }
  
  if (!is.numeric(rel_change_tol) || length(rel_change_tol) != 1 || 
      rel_change_tol <= 0) {
    stop("rel_change_tol must be a positive scalar")
  }
  
  # Validate each condition design matrix
  for (c in 1:k) {
    if (!is.matrix(X_condition_list_proj_matrices[[c]])) {
      stop(sprintf("X_condition_list_proj_matrices[[%d]] must be a matrix", c))
    }
    if (nrow(X_condition_list_proj_matrices[[c]]) != n) {
      stop(sprintf(
        "X_condition_list_proj_matrices[[%d]] must have %d rows to match Y_proj_matrix", 
        c, n
      ))
    }
    if (ncol(X_condition_list_proj_matrices[[c]]) != p) {
      stop(sprintf(
        "X_condition_list_proj_matrices[[%d]] must have %d columns to match H_shapes_allvox_matrix rows", 
        c, p
      ))
    }
  }
  
  # Initialize output matrix
  Beta_condition_final_matrix <- matrix(0, nrow = k, ncol = V)
  
  # For MVP, we only do one iteration (max_iter = 1)
  # Future versions could iterate between beta and HRF estimation
  for (iter in 1:max_iter) {
    
    # Store previous beta for convergence check (if max_iter > 1)
    if (iter > 1) {
      Beta_previous <- Beta_condition_final_matrix
    }
    
    # Main voxel loop
    # This could be parallelized in future versions
    for (v in 1:V) {
      
      # Extract data and HRF for current voxel
      Y_voxel <- Y_proj_matrix[, v]
      H_voxel <- H_shapes_allvox_matrix[, v]
      
      # Construct design matrix for this voxel
      # Each column is a condition-specific regressor convolved with voxel HRF
      X_design_voxel <- matrix(0, nrow = n, ncol = k)
      
      for (c in 1:k) {
        # Convolve condition design with voxel-specific HRF
        # X_c is n x p, H_voxel is p x 1, result is n x 1
        X_design_voxel[, c] <- X_condition_list_proj_matrices[[c]] %*% H_voxel
      }
      
      # Ridge regression for this voxel
      # beta = (X'X + lambda*I)^(-1) * X'y
      XtX_voxel <- crossprod(X_design_voxel)  # k x k
      XtX_voxel_reg <- XtX_voxel + lambda_beta_final * diag(k)
      XtY_voxel <- crossprod(X_design_voxel, Y_voxel)  # k x 1
      
      # Solve for beta
      # Check for numerical issues
      if (any(is.na(XtX_voxel_reg)) || any(is.infinite(XtX_voxel_reg))) {
        warning(sprintf("Numerical issues in voxel %d, setting betas to zero", v))
        Beta_condition_final_matrix[, v] <- 0
      } else {
        tryCatch({
          beta_voxel <- solve(XtX_voxel_reg, XtY_voxel)
          Beta_condition_final_matrix[, v] <- as.vector(beta_voxel)
        }, error = function(e) {
          warning(sprintf("Failed to solve for voxel %d: %s. Setting betas to zero.", 
                         v, e$message))
          Beta_condition_final_matrix[, v] <- 0
        })
      }
    }
    
    # Check convergence (only relevant if max_iter > 1)
    if (iter > 1) {
      # Compute relative change in beta estimates
      beta_change <- norm(Beta_condition_final_matrix - Beta_previous, "F") / 
                     max(norm(Beta_previous, "F"), .Machine$double.eps)
      
      if (beta_change < rel_change_tol) {
        message(sprintf("Converged after %d iterations (relative change = %.6f)", 
                       iter, beta_change))
        break
      }
    }
    
    # If max_iter > 1, we could re-estimate HRFs here using current betas
    # This would implement the full alternating optimization
    # For MVP (max_iter = 1), we skip this
  }
  
  return(Beta_condition_final_matrix)
}