# Improved Core Voxel-wise Fit Engine Functions (Component 1)
# Addresses numerical stability, performance, and architectural issues

#' VoxelWiseGLM R6 Class for Efficient Voxel-wise Analysis
#'
#' An object-oriented implementation of the voxel-wise GLM fitting engine
#' with pre-computation of shared matrices and support for parallelization.
#'
#' @field Z Confound matrix (n x q)
#' @field Z_qr QR decomposition of confound matrix Z
#' @field Q_Z Q matrix from QR decomposition of Z
#' @field rank_Z Rank of the confound matrix Z
#' @field n_timepoints Number of timepoints in the data
#' @field ridge_lambda Ridge regularization parameter (default: 1e-6)
#' @field use_parallel Whether to use parallel processing (default: TRUE)
#' @field chunk_size Number of voxels to process per chunk (default: 1000)
#'
#' @importFrom R6 R6Class
#' @importFrom Matrix sparseMatrix Diagonal
#' @export
VoxelWiseGLM <- R6::R6Class(
  "VoxelWiseGLM",
  public = list(
    # Pre-computed matrices
    Z = NULL,           # Confound matrix
    Z_qr = NULL,        # QR decomposition of Z
    Q_Z = NULL,         # Q from QR(Z)
    rank_Z = NULL,      # Rank of Z
    n_timepoints = NULL,
    
    # Configuration
    ridge_lambda = 1e-6,
    use_parallel = TRUE,
    chunk_size = 1000,
    
    #' Initialize the GLM engine
    #' @param confounds Optional n x q confound matrix
    #' @param ridge_lambda Ridge regularization parameter
    initialize = function(confounds = NULL, ridge_lambda = 1e-6) {
      self$ridge_lambda <- ridge_lambda
      
      if (!is.null(confounds)) {
        if (!is.matrix(confounds)) {
          stop("confounds must be a matrix or NULL")
        }
        
        self$Z <- confounds
        self$n_timepoints <- nrow(confounds)
        
        # QR decomposition with column pivoting for rank detection
        self$Z_qr <- qr(self$Z, LAPACK = TRUE)
        self$rank_Z <- self$Z_qr$rank
        
        if (self$rank_Z < ncol(self$Z)) {
          warning(sprintf(
            "Confound matrix is rank deficient: rank %d < %d columns. Using %d independent components.",
            self$rank_Z, ncol(self$Z), self$rank_Z
          ))
        }
        
        # Extract Q matrix (only rank columns)
        self$Q_Z <- qr.Q(self$Z_qr)[, seq_len(self$rank_Z), drop = FALSE]
      }
    },
    
    #' Project out confounds from data and design matrices
    #' @param Y n x V data matrix
    #' @param X_list List of n x p design matrices
    #' @return List with projected Y and X matrices
    project_out_confounds = function(Y, X_list = NULL) {
      if (is.null(self$Q_Z)) {
        # No confounds to project
        return(list(Y_proj = Y, X_proj = X_list))
      }
      
      # Efficient projection: Y_proj = Y - Q(Q'Y)
      # Using crossprod for efficiency
      QtY <- crossprod(self$Q_Z, Y)
      Y_proj <- Y - self$Q_Z %*% QtY
      
      X_proj <- NULL
      if (!is.null(X_list)) {
        X_proj <- lapply(X_list, function(X) {
          QtX <- crossprod(self$Q_Z, X)
          X - self$Q_Z %*% QtX
        })
      }
      
      return(list(Y_proj = Y_proj, X_proj = X_proj))
    },
    
    #' Fit GLM with QR decomposition
    #' @param Y n x V data matrix
    #' @param X n x k design matrix or list of matrices
    #' @param project_confounds Whether to project out confounds first
    #' @return Coefficient matrix (k x V)
    fit = function(Y, X, project_confounds = TRUE) {
      # Handle list input
      if (is.list(X)) {
        X <- do.call(cbind, X)
      }
      
      # Validate dimensions
      n <- nrow(Y)
      V <- ncol(Y)
      k <- ncol(X)
      
      if (nrow(X) != n) {
        stop("X and Y must have same number of rows")
      }
      
      # Project out confounds if requested
      if (project_confounds && !is.null(self$Q_Z)) {
        proj_data <- self$project_out_confounds(Y, list(X))
        Y <- proj_data$Y_proj
        X <- proj_data$X_proj[[1]]
      }
      
      # Fit using QR decomposition (more stable than normal equations)
      if (self$ridge_lambda > 0) {
        # Ridge regression via augmented system
        X_aug <- rbind(X, sqrt(self$ridge_lambda) * diag(k))
        Y_aug <- rbind(Y, matrix(0, k, V))
        coef <- qr.solve(X_aug, Y_aug)
      } else {
        # Standard least squares via QR
        coef <- qr.solve(X, Y)
      }
      
      return(coef)
    },
    
    #' Fit with parallelization over voxels
    #' @param Y n x V data matrix
    #' @param X n x k design matrix
    #' @param n_cores Number of cores to use
    #' @return Coefficient matrix
    fit_parallel = function(Y, X, n_cores = NULL) {
      if (!requireNamespace("future.apply", quietly = TRUE)) {
        warning("future.apply not available, using sequential fit")
        return(self$fit(Y, X))
      }
      
      V <- ncol(Y)
      
      # Determine chunk assignments
      n_chunks <- ceiling(V / self$chunk_size)
      chunk_indices <- split(seq_len(V), 
                           rep(seq_len(n_chunks), 
                               length.out = V))
      
      # Process chunks in parallel
      results <- future.apply::future_lapply(chunk_indices, function(idx) {
        self$fit(Y[, idx, drop = FALSE], X)
      })
      
      # Combine results
      do.call(cbind, results)
    }
  )
)

#' Transform designs to manifold basis (improved)
#'
#' Vectorized transformation of design matrices to manifold basis
#'
#' @param X_list List of n x p design matrices  
#' @param B p x m manifold basis matrix
#' @return List of n x m transformed design matrices
#' @export
transform_designs_to_manifold_basis_improved <- function(X_list, B) {
  if (!is.list(X_list)) {
    stop("X_list must be a list of matrices")
  }
  
  if (!is.matrix(B)) {
    stop("B must be a matrix")
  }
  
  # Validate dimensions
  p <- nrow(B)
  
  # Transform each design matrix
  lapply(X_list, function(X) {
    if (ncol(X) != p) {
      stop(sprintf("Design matrix has %d columns but B has %d rows", 
                   ncol(X), p))
    }
    X %*% B
  })
}

#' Extract manifold coordinates and amplitudes via SVD (vectorized)
#'
#' Efficient extraction using block processing to avoid R loops
#'
#' @param Gamma (k*m) x V coefficient matrix
#' @param m Manifold dimension
#' @param k Number of conditions
#' @param block_size Number of voxels to process at once
#' @return List with Xi (m x V) and Beta (k x V) matrices
#' @export
extract_xi_beta_svd_block <- function(Gamma, m, k, block_size = 100) {
  V <- ncol(Gamma)
  
  if (nrow(Gamma) != k * m) {
    stop("Gamma dimensions inconsistent with k and m")
  }
  
  # Pre-allocate output
  Xi <- matrix(0, m, V)
  Beta <- matrix(0, k, V)
  
  # Process in blocks for efficiency
  for (start in seq(1, V, by = block_size)) {
    end <- min(start + block_size - 1, V)
    idx <- start:end
    
    # Extract block of gamma coefficients
    gamma_block <- Gamma[, idx, drop = FALSE]
    
    # Process each voxel in block
    for (i in seq_along(idx)) {
      v <- idx[i]
      
      # Reshape to m x k matrix (consistent orientation)
      G_v <- matrix(gamma_block[, i], nrow = m, ncol = k)
      
      # Skip if all zeros
      if (all(abs(G_v) < .Machine$double.eps)) {
        next
      }
      
      # Compute SVD
      svd_result <- tryCatch({
        svd(G_v, nu = 1, nv = 1)
      }, error = function(e) {
        warning(sprintf("SVD failed for voxel %d: %s", v, e$message))
        NULL
      })
      
      if (!is.null(svd_result) && length(svd_result$d) > 0) {
        # Extract first singular vectors scaled by singular value
        Xi[, v] <- svd_result$u[, 1] * sqrt(svd_result$d[1])
        Beta[, v] <- svd_result$v[, 1] * sqrt(svd_result$d[1])
      }
    }
  }
  
  return(list(Xi_raw_matrix = Xi, Beta_raw_matrix = Beta))
}

#' Apply intrinsic identifiability constraints (vectorized)
#'
#' Efficient sign/scale alignment without redundant computations
#'
#' @param Xi_raw m x V raw manifold coordinates
#' @param Beta_raw k x V raw amplitudes  
#' @param B p x m manifold basis
#' @param h_ref p x 1 reference HRF shape
#' @param h_mode Sign alignment mode
#' @param scale_method Method for scaling: "l2_norm", "max_abs_val", "beta_norm", or "none"
#' @return List with aligned Xi and Beta
#' @export
apply_identifiability_vectorized <- function(Xi_raw, Beta_raw, B, h_ref, 
                                           h_mode = "max_correlation",
                                           scale_method = "l2_norm") {
  m <- nrow(Xi_raw)
  k <- nrow(Beta_raw)
  V <- ncol(Xi_raw)
  
  # Pre-compute reference manifold coordinates once
  xi_ref <- as.vector(MASS::ginv(B) %*% h_ref)
  
  # Pre-allocate output
  Xi_ident <- Xi_raw
  Beta_ident <- Beta_raw
  
  # Vectorized operations where possible
  if (h_mode == "max_abs") {
    # Find dimension with max absolute reference value
    max_dim <- which.max(abs(xi_ref))
    
    # Vectorized sign flip
    signs <- sign(Xi_raw[max_dim, ])
    signs[signs == 0] <- 1
    
    Xi_ident <- Xi_raw * rep(signs, each = m)
    Beta_ident <- Beta_raw * rep(signs, each = k)
    
  } else if (h_mode == "max_correlation") {
    # Compute correlations efficiently
    xi_ref_norm <- xi_ref / sqrt(sum(xi_ref^2))
    
    # Normalize Xi columns
    Xi_norms <- sqrt(colSums(Xi_raw^2))
    Xi_norm <- Xi_raw / rep(Xi_norms, each = m)
    Xi_norm[, Xi_norms == 0] <- 0
    
    # Compute correlations via matrix multiplication
    correlations <- as.vector(crossprod(xi_ref_norm, Xi_norm))
    
    # Apply sign based on correlation
    signs <- sign(correlations)
    signs[signs == 0] <- 1
    
    Xi_ident <- Xi_raw * rep(signs, each = m)
    Beta_ident <- Beta_raw * rep(signs, each = k)
    
  } else if (h_mode == "first_positive") {
    # Find first non-zero element in each column
    for (v in 1:V) {
      first_nonzero <- which(abs(Xi_raw[, v]) > 1e-10)[1]
      if (!is.na(first_nonzero) && Xi_raw[first_nonzero, v] < 0) {
        Xi_ident[, v] <- -Xi_raw[, v]
        Beta_ident[, v] <- -Beta_raw[, v]
      }
    }
  }
  
  # Apply scaling based on method
  if (scale_method == "l2_norm") {
    # L2 normalization of HRF 
    # To preserve signal (Xi * Beta), we need to redistribute the scaling
    for (v in 1:V) {
      # Check if Xi values are effectively zero (tiny values)
      if (all(abs(Xi_ident[, v]) < 1e-10)) {
        Xi_ident[, v] <- 0
        Beta_ident[, v] <- 0
      } else if (any(Xi_ident[, v] != 0)) {
        hrf <- B %*% Xi_ident[, v]
        hrf_norm <- sqrt(sum(hrf^2))
        if (hrf_norm > 0) {
          # Scale Xi to normalize HRF
          Xi_ident[, v] <- Xi_ident[, v] / hrf_norm
          # Scale Beta inversely to preserve signal
          Beta_ident[, v] <- Beta_ident[, v] * hrf_norm
        }
      }
    }
  } else if (scale_method == "max_abs_val") {
    # Max absolute value normalization of HRF
    for (v in 1:V) {
      # Check if Xi values are effectively zero (tiny values)
      if (all(abs(Xi_ident[, v]) < 1e-10)) {
        Xi_ident[, v] <- 0
        Beta_ident[, v] <- 0
      } else if (any(Xi_ident[, v] != 0)) {
        hrf <- B %*% Xi_ident[, v]
        max_abs <- max(abs(hrf))
        if (max_abs > 0) {
          # Scale Xi to normalize HRF
          Xi_ident[, v] <- Xi_ident[, v] / max_abs
          # Scale Beta inversely to preserve signal
          Beta_ident[, v] <- Beta_ident[, v] * max_abs
        }
      }
    }
  } else if (scale_method == "beta_norm") {
    # Original Beta normalization (for backward compatibility)
    # First check for tiny Xi values
    for (v in 1:V) {
      if (all(abs(Xi_ident[, v]) < 1e-10)) {
        Xi_ident[, v] <- 0
        Beta_ident[, v] <- 0
      }
    }
    
    beta_norms <- sqrt(colSums(Beta_ident^2))
    beta_norms[beta_norms == 0] <- 1
    
    Xi_ident <- Xi_ident * rep(beta_norms, each = m)
    Beta_ident <- Beta_ident / rep(beta_norms, each = k)
  } else {
    # scale_method == "none", but still check for tiny values
    for (v in 1:V) {
      if (all(abs(Xi_ident[, v]) < 1e-10)) {
        Xi_ident[, v] <- 0
        Beta_ident[, v] <- 0
      }
    }
  }
  
  return(list(
    Xi_ident_matrix = Xi_ident,
    Beta_ident_matrix = Beta_ident
  ))
}

# Wrapper functions for backward compatibility

#' Project out confounds (backward compatible wrapper)
#' @export
project_out_confounds_core <- function(Y_data_matrix, X_list_of_matrices, 
                                      Z_confounds_matrix = NULL) {
  engine <- VoxelWiseGLM$new(confounds = Z_confounds_matrix)
  engine$project_out_confounds(Y_data_matrix, X_list_of_matrices)
}

#' Solve GLM for gamma (improved wrapper)
#' @export
solve_glm_for_gamma_core <- function(Z_list_of_matrices, Y_proj_matrix,
                                    lambda_gamma = 0, orthogonal_approx_flag = FALSE) {
  if (orthogonal_approx_flag) {
    warning("orthogonal_approx_flag is deprecated and ignored")
  }
  
  engine <- VoxelWiseGLM$new(ridge_lambda = lambda_gamma)
  X <- do.call(cbind, Z_list_of_matrices)
  engine$fit(Y_proj_matrix, X, project_confounds = FALSE)
}

#' Extract SVD components (wrapper)
#' @export
extract_xi_beta_raw_svd_core <- function(Gamma_coeffs_matrix, m_manifold_dim, 
                                        k_conditions) {
  extract_xi_beta_svd_block(Gamma_coeffs_matrix, m_manifold_dim, k_conditions)
}