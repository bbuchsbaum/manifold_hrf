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
#' @param n_jobs Number of parallel workers for voxel-wise computation.
#'   
#' @return Beta_condition_final_matrix A k x V matrix of final condition-level 
#'   beta estimates
#'   
#' @details This function implements Component 4 of the M-HRF-LSS pipeline.
#'   It re-estimates condition-level amplitudes using the final smoothed HRF
#'   shapes. 
#'   
#'   Memory optimization: Instead of precomputing all k × (n × V) matrices
#'   (which requires ~10GB for whole brain), this implementation computes
#'   regressors on-the-fly for each voxel, reducing memory usage to O(n × k).
#'   
#'   Numerical stability: Uses Cholesky decomposition when possible, falling
#'   back to QR decomposition for non-positive definite cases.
#'   
#'   The MVP version does a single pass (max_iter=1), but the framework 
#'   supports iterative refinement between HRFs and betas.
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
                                                                    rel_change_tol = 1e-4),
                                              n_jobs = 1) {
  
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
  
  lambda_beta_final <- .validate_and_standardize_lambda(lambda_beta_final,
                                                        "lambda_beta_final")
  
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
    
    voxel_fun <- function(v) {
      # Build n x k design matrix on the fly for this voxel
      D_v <- vapply(X_condition_list_proj_matrices,
                    function(Xc) Xc %*% H_shapes_allvox_matrix[, v],
                    numeric(n))
      
      # Compute cross products for this voxel
      XtX <- crossprod(D_v) + lambda_beta_final * diag(k)
      XtY <- crossprod(D_v, Y_proj_matrix[, v])
      
      # Use Cholesky decomposition for stability
      tryCatch({
        chol_XtX <- chol(XtX)
        betas <- backsolve(chol_XtX, forwardsolve(t(chol_XtX), XtY))
        as.vector(betas)
      }, error = function(e) {
        # Fall back to QR if Cholesky fails (matrix not positive definite)
        tryCatch({
          qr_D <- qr(D_v)
          betas <- qr.coef(qr_D, Y_proj_matrix[, v])
          as.vector(betas)
        }, error = function(e2) {
          warning(sprintf("Failed to solve for voxel %d: %s. Setting betas to zero.",
                          v, e2$message))
          rep(0, k)
        })
      })
    }

    # Use centralized parallel processing helper if available
    if (exists(".lss_process_voxels", mode = "function")) {
      res_list <- .lss_process_voxels(
        voxel_indices = seq_len(V),
        worker_fun = voxel_fun,
        n_cores = n_jobs,
        progress = FALSE,
        .globals = c("X_condition_list_proj_matrices", "H_shapes_allvox_matrix", 
                     "Y_proj_matrix", "lambda_beta_final", "n", "k")
      )
    } else {
      # Fall back to .parallel_lapply if helper not available
      res_list <- .parallel_lapply(seq_len(V), voxel_fun, n_jobs)
    }
    Beta_condition_final_matrix <- do.call(cbind, res_list)
    
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