# Robust Voxelwise Fitting Improvements
# Implementation of SOUND-VOXFIT-REGULARIZE

#' Robust SVD Extraction with Conditioning
#'
#' Enhanced version of extract_xi_beta_raw_svd_core with automatic regularization
#' and fallback strategies
#'
#' @param Gamma_coeffs_matrix (km) x V matrix of stacked gamma coefficients
#' @param m_manifold_dim Number of manifold dimensions
#' @param k_conditions Number of conditions
#' @param regularization_factor Factor to increase regularization on poor conditioning
#' @param max_condition_number Maximum acceptable condition number
#' @param use_randomized_svd Use randomized SVD for large problems
#' @return List with Xi_raw_matrix, Beta_raw_matrix, and quality metrics
#' @export
extract_xi_beta_raw_svd_robust <- function(Gamma_coeffs_matrix,
                                          m_manifold_dim,
                                          k_conditions,
                                          regularization_factor = 10,
                                          max_condition_number = 1e8,
                                          use_randomized_svd = FALSE,
                                          logger = NULL) {
  
  # Get dimensions
  km <- nrow(Gamma_coeffs_matrix)
  V <- ncol(Gamma_coeffs_matrix)
  
  if (km != k_conditions * m_manifold_dim) {
    stop("Gamma_coeffs_matrix has incorrect number of rows")
  }
  
  # Initialize output
  Xi_raw <- matrix(0, m_manifold_dim, V)
  Beta_raw <- matrix(0, k_conditions, V)
  
  # Track quality metrics
  quality_metrics <- list(
    condition_numbers = numeric(V),
    svd_method = character(V),
    regularization_applied = logical(V),
    singular_value_gaps = numeric(V)
  )
  
  # Process each voxel
  for (v in 1:V) {
    
    # Reshape gamma for this voxel
    gamma_v <- Gamma_coeffs_matrix[, v]
    Gamma_mat <- matrix(gamma_v, nrow = k_conditions, ncol = m_manifold_dim, byrow = TRUE)
    
    # Check if gamma is all zeros
    if (all(abs(gamma_v) < .Machine$double.eps)) {
      Xi_raw[, v] <- 0
      Beta_raw[, v] <- 0
      quality_metrics$svd_method[v] <- "zero"
      next
    }
    
    # Compute condition number
    gamma_scale <- max(abs(Gamma_mat))
    if (gamma_scale > 0) {
      Gamma_scaled <- Gamma_mat / gamma_scale
      cn <- kappa(Gamma_scaled)
      quality_metrics$condition_numbers[v] <- cn
      
      # Check if matrix is rank-1 (special case)
      svd_check <- svd(Gamma_mat, nu = 0, nv = 0)
      n_sig <- sum(svd_check$d > max(svd_check$d) * .Machine$double.eps * 100)
      
      # If poorly conditioned but NOT rank-1, add regularization
      if (cn > max_condition_number && n_sig > 1) {
        reg_amount <- (cn / max_condition_number) * regularization_factor * .Machine$double.eps
        # Add regularization to diagonal
        diag(Gamma_mat) <- diag(Gamma_mat) + reg_amount
        quality_metrics$regularization_applied[v] <- TRUE
        warning(sprintf("Voxel %d: Applied regularization due to condition number %.2e", v, cn))
      }
    }
    
    # Try SVD with error handling
    svd_result <- tryCatch({
      
      if (use_randomized_svd && k_conditions > 10 && m_manifold_dim > 10) {
        # Use randomized SVD for efficiency
        quality_metrics$svd_method[v] <- "randomized"
        if (requireNamespace("rsvd", quietly = TRUE)) {
          rsvd::rsvd(Gamma_mat, k = min(k_conditions, m_manifold_dim))
        } else {
          svd(Gamma_mat)
        }
      } else {
        quality_metrics$svd_method[v] <- "standard"
        svd(Gamma_mat)
      }
      
    }, error = function(e) {
      warning(sprintf("SVD failed for voxel %d: %s. Using fallback.", v, e$message))
      quality_metrics$svd_method[v] <- "fallback"
      # Fallback: use first PC or zeros if degenerate
      if (all(is.finite(Gamma_mat)) && sum(Gamma_mat^2) > 0) {
        list(
          u = matrix(1/sqrt(k_conditions), k_conditions, 1),
          v = matrix(1/sqrt(m_manifold_dim), m_manifold_dim, 1),
          d = sqrt(sum(Gamma_mat^2))
        )
      } else {
        # Complete fallback for degenerate case
        list(
          u = matrix(0, k_conditions, 1),
          v = matrix(0, m_manifold_dim, 1),
          d = 0
        )
      }
    })
    
    # Extract components with "soft" truncation based on singular value gap
    d <- svd_result$d
    
    if (length(d) > 1) {
      # Find significant gap in singular values
      gaps <- diff(d) / d[-length(d)]
      quality_metrics$singular_value_gaps[v] <- max(abs(gaps))
      
      # Soft truncation: weight by relative singular value
      weights <- d / (d[1] + .Machine$double.eps)
      weights[weights < 0.01] <- 0  # Threshold at 1% of largest
    } else {
      weights <- 1
    }
    
    # Extract xi and beta
    if (length(d) > 0 && d[1] > .Machine$double.eps) {
      # Weighted reconstruction
      Xi_raw[, v] <- svd_result$v[, 1] * sqrt(d[1]) * weights[1]
      Beta_raw[, v] <- svd_result$u[, 1] * sqrt(d[1]) * weights[1]
    } else {
      # Degenerate case
      Xi_raw[, v] <- 0
      Beta_raw[, v] <- 0
      quality_metrics$svd_method[v] <- "degenerate"
    }
  }
  
  # Report quality summary
  n_regularized <- sum(quality_metrics$regularization_applied)
  if (n_regularized > 0) {
    msg <- sprintf("Applied regularization to %d/%d voxels", n_regularized, V)
    if (!is.null(logger)) logger$add(msg) else message(msg)
  }

  n_fallback <- sum(quality_metrics$svd_method == "fallback")
  if (n_fallback > 0) {
    msg <- sprintf("Used fallback SVD for %d/%d voxels", n_fallback, V)
    if (!is.null(logger)) logger$add(msg) else message(msg)
  }
  
  return(list(
    Xi_raw_matrix = Xi_raw,
    Beta_raw_matrix = Beta_raw,
    quality_metrics = quality_metrics
  ))
}


#' Smart Initialization for Voxelwise Fit
#'
#' Provides intelligent starting values for the manifold fit
#'
#' @param Y_data n x V data matrix
#' @param X_condition_list List of condition design matrices
#' @param hrf_canonical Canonical HRF for initialization
#' @param use_spatial_clusters Whether to use spatial clustering
#' @param voxel_coords Voxel coordinates for clustering
#' @param m_manifold_dim Manifold dimension for Xi initialization
#' @return List with initial Xi and Beta matrices
#' @export
smart_initialize <- function(Y_data, X_condition_list, hrf_canonical,
                           use_spatial_clusters = TRUE,
                           voxel_coords = NULL,
                           m_manifold_dim = 5) {
  
  n <- nrow(Y_data)
  V <- ncol(Y_data)
  k <- length(X_condition_list)
  
  message("Computing smart initialization...")
  
  # First do standard GLM with canonical HRF
  X_canonical <- matrix(0, n, k)
  for (c in 1:k) {
    X_canonical[, c] <- X_condition_list[[c]] %*% hrf_canonical
  }
  
  # Solve for initial betas
  XtX <- crossprod(X_canonical)
  XtX_reg <- XtX + diag(0.01, k)  # Small regularization
  XtX_inv <- solve(XtX_reg)
  
  Beta_init <- matrix(0, k, V)
  R2_init <- numeric(V)
  
  for (v in 1:V) {
    beta_v <- XtX_inv %*% crossprod(X_canonical, Y_data[, v])
    Beta_init[, v] <- beta_v
    
    # Compute R² for quality assessment
    y_pred <- X_canonical %*% beta_v
    ss_tot <- sum((Y_data[, v] - mean(Y_data[, v]))^2)
    ss_res <- sum((Y_data[, v] - y_pred)^2)
    
    # Handle edge cases: if ss_tot is 0 (constant data), R² is undefined
    if (ss_tot < .Machine$double.eps) {
      R2_init[v] <- 0  # Constant data has no variance to explain
    } else {
      R2_init[v] <- 1 - ss_res / ss_tot
      # Clamp to [0, 1] to handle numerical issues
      R2_init[v] <- max(0, min(1, R2_init[v]))
    }
  }
  
  # Identify well-fit voxels to use as exemplars
  good_voxels <- which(R2_init > quantile(R2_init, 0.75, na.rm = TRUE))
  
  # If using spatial clustering, find similar voxels
  if (use_spatial_clusters && !is.null(voxel_coords) && length(good_voxels) > 10) {
    
    # Cluster good voxels
    n_clusters <- min(20, length(good_voxels) / 5)
    kmeans_result <- kmeans(voxel_coords[good_voxels, ], centers = n_clusters)
    
    # For each voxel, find nearest cluster center
    cluster_centers <- kmeans_result$centers
    nearest_cluster <- apply(voxel_coords, 1, function(v) {
      distances <- rowSums((cluster_centers - matrix(v, n_clusters, 3, byrow = TRUE))^2)
      which.min(distances)
    })
    
    message(sprintf("Using %d spatial clusters for initialization", n_clusters))
  } else {
    # No spatial clustering - use global exemplar
    nearest_cluster <- rep(1, V)
  }
  
  # Create initial Xi as small perturbations
  # This helps avoid local minima
  Xi_init <- matrix(rnorm(m_manifold_dim * V, sd = 0.1), m_manifold_dim, V)
  
  # Scale Xi based on initial beta magnitudes
  beta_scale <- apply(Beta_init, 2, function(x) sqrt(sum(x^2)))
  Xi_init <- Xi_init * rep(beta_scale, each = m_manifold_dim)
  
  return(list(
    Xi_init = Xi_init,
    Beta_init = Beta_init,
    R2_init = R2_init,
    good_voxels = good_voxels,
    nearest_cluster = nearest_cluster
  ))
}


#' Check and Fix Stuck Voxels
#'
#' Detects voxels with very low variance in Xi and reinitializes them
#'
#' @param Xi_current Current Xi estimates
#' @param Xi_previous Previous Xi estimates
#' @param variance_threshold Threshold for detecting stuck voxels
#' @param reinit_sd Standard deviation for reinitialization
#' @return Updated Xi matrix
#' @keywords internal
fix_stuck_voxels <- function(Xi_current, Xi_previous, 
                           variance_threshold = 1e-6,
                           reinit_sd = 0.1) {
  
  m <- nrow(Xi_current)
  V <- ncol(Xi_current)
  
  # Compute change from previous iteration
  if (!is.null(Xi_previous)) {
    Xi_change <- Xi_current - Xi_previous
    change_variance <- apply(Xi_change, 2, var)
    
    # Identify stuck voxels
    stuck <- which(change_variance < variance_threshold)
    
    if (length(stuck) > 0) {
      message(sprintf("Reinitializing %d stuck voxels", length(stuck)))
      
      # Reinitialize with small random values
      Xi_current[, stuck] <- matrix(rnorm(m * length(stuck), sd = reinit_sd),
                                   m, length(stuck))
    }
  }
  
  return(Xi_current)
}
