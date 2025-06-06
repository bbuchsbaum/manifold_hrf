# Robust Spatial Smoothing and Outlier Handling
# Implementation of SOUND-SPATIAL-ADAPTIVE and SOUND-OUTLIER-ROBUST

#' Compute Local SNR for Adaptive Smoothing
#'
#' Estimates signal-to-noise ratio for each voxel to guide smoothing strength
#'
#' @param Y_data n x V data matrix
#' @param Y_predicted n x V predicted data (optional)
#' @param method Method for SNR estimation
#' @return Vector of SNR values for each voxel
#' @keywords internal
compute_local_snr <- function(Y_data, Y_predicted = NULL, 
                             method = c("temporal_variance", "residual")) {
  
  method <- match.arg(method)
  V <- ncol(Y_data)
  snr <- numeric(V)
  
  if (method == "temporal_variance") {
    # SNR based on temporal variance
    for (v in 1:V) {
      y <- Y_data[, v]
      # Robust estimate using MAD
      signal_var <- var(y)
      noise_mad <- median(abs(diff(y))) * 1.4826
      noise_var <- noise_mad^2
      
      snr[v] <- signal_var / (noise_var + .Machine$double.eps)
    }
  } else if (method == "residual" && !is.null(Y_predicted)) {
    # SNR based on fit residuals
    for (v in 1:V) {
      signal_var <- var(Y_predicted[, v])
      residuals <- Y_data[, v] - Y_predicted[, v]
      noise_var <- var(residuals)
      
      snr[v] <- signal_var / (noise_var + .Machine$double.eps)
    }
  }
  
  # Cap extreme values
  snr[is.na(snr)] <- 1
  snr[is.infinite(snr)] <- 100
  snr[snr > 100] <- 100
  snr[snr < 0.1] <- 0.1
  
  return(snr)
}


#' Adaptive Spatial Smoothing
#'
#' Enhanced spatial smoothing with SNR-adaptive regularization
#'
#' @param Xi_ident_matrix m x V manifold coordinates
#' @param L_sp_sparse_matrix V x V spatial Laplacian
#' @param lambda_spatial_smooth Base smoothing parameter
#' @param local_snr Vector of SNR values for adaptive smoothing
#' @param edge_preserve Whether to preserve edges
#' @param voxel_coords Voxel coordinates for edge detection
#' @return Smoothed Xi matrix
#' @export
apply_spatial_smoothing_adaptive <- function(Xi_ident_matrix,
                                           L_sp_sparse_matrix,
                                           lambda_spatial_smooth,
                                           local_snr = NULL,
                                           edge_preserve = TRUE,
                                           voxel_coords = NULL) {
  
  m <- nrow(Xi_ident_matrix)
  V <- ncol(Xi_ident_matrix)
  
  # If no SNR provided, use uniform smoothing
  if (is.null(local_snr)) {
    local_snr <- rep(1, V)
  }
  
  # Adaptive lambda based on SNR
  # Low SNR → more smoothing (higher lambda)
  # High SNR → less smoothing (lower lambda)
  snr_factor <- 1 / sqrt(local_snr)
  snr_factor <- snr_factor / median(snr_factor)  # Normalize
  
  lambda_adaptive <- lambda_spatial_smooth * snr_factor
  
  # Edge detection if requested
  if (edge_preserve && !is.null(voxel_coords)) {
    edge_weights <- compute_edge_weights(Xi_ident_matrix, voxel_coords)
  } else {
    edge_weights <- rep(1, V)
  }
  
  # Apply adaptive smoothing
  Xi_smoothed <- matrix(0, m, V)
  
  # Check for isolated voxels
  if (inherits(L_sp_sparse_matrix, "Matrix")) {
    row_sums <- Matrix::rowSums(abs(L_sp_sparse_matrix))
  } else {
    row_sums <- rowSums(abs(L_sp_sparse_matrix))
  }
  isolated_voxels <- which(row_sums < .Machine$double.eps)
  
  if (length(isolated_voxels) > 0) {
    message(sprintf("Skipping smoothing for %d isolated voxels", length(isolated_voxels)))
  }
  
  # Smooth each manifold dimension
  for (j in 1:m) {
    xi_j <- Xi_ident_matrix[j, ]
    
    # Create adaptive system matrix
    # (I + λ_adaptive * L * edge_weights)
    if (inherits(L_sp_sparse_matrix, "Matrix")) {
      # Sparse implementation
      Lambda_diag <- Matrix::Diagonal(x = lambda_adaptive * edge_weights)
      A_matrix <- Matrix::Diagonal(V) + Lambda_diag %*% L_sp_sparse_matrix
      
      # Solve with error handling
      xi_smooth <- tryCatch({
        Matrix::solve(A_matrix, xi_j)
      }, error = function(e) {
        warning(sprintf("Adaptive smoothing failed for dimension %d: %s", j, e$message))
        xi_j  # Return original
      })
    } else {
      # Dense implementation
      A_matrix <- diag(V) + diag(lambda_adaptive * edge_weights) %*% L_sp_sparse_matrix
      
      xi_smooth <- tryCatch({
        solve(A_matrix, xi_j)
      }, error = function(e) {
        warning(sprintf("Adaptive smoothing failed for dimension %d: %s", j, e$message))
        xi_j  # Return original
      })
    }
    
    # Keep isolated voxels unchanged
    xi_smooth[isolated_voxels] <- xi_j[isolated_voxels]
    
    Xi_smoothed[j, ] <- xi_smooth
  }
  
  # Report adaptive smoothing summary
  message(sprintf("Adaptive smoothing applied: lambda range [%.3f, %.3f]",
                 min(lambda_adaptive), max(lambda_adaptive)))
  
  return(Xi_smoothed)
}


#' Compute Edge Weights for Edge-Preserving Smoothing
#'
#' Detects edges in manifold coordinates to reduce smoothing across boundaries
#'
#' @param Xi_matrix m x V manifold coordinates
#' @param voxel_coords V x 3 spatial coordinates
#' @param edge_threshold Threshold for edge detection
#' @return Vector of edge weights (0 = edge, 1 = smooth region)
#' @keywords internal
compute_edge_weights <- function(Xi_matrix, voxel_coords, edge_threshold = 2) {
  
  V <- ncol(Xi_matrix)
  edge_weights <- rep(1, V)
  
  # Compute local gradient magnitude
  for (v in 1:V) {
    # Find spatial neighbors
    # Fix array recycling: use rep() to properly expand the vector
    voxel_v_matrix <- matrix(rep(voxel_coords[v, ], each = V), nrow = V, ncol = 3)
    distances <- sqrt(rowSums((voxel_coords - voxel_v_matrix)^2))
    neighbors <- which(distances > 0 & distances < 2)
    
    if (length(neighbors) > 0) {
      # Compute gradient in manifold space
      xi_v <- Xi_matrix[, v]
      xi_neighbors <- Xi_matrix[, neighbors, drop = FALSE]
      
      # Mean squared difference to neighbors
      gradient_mag <- mean(colSums((xi_neighbors - xi_v)^2))
      
      # Convert to weight (high gradient → low weight)
      if (gradient_mag > edge_threshold) {
        edge_weights[v] <- exp(-gradient_mag / edge_threshold)
      }
    }
  }
  
  return(edge_weights)
}


#' Detect and Handle Outlier Timepoints
#'
#' Identifies outlier timepoints and returns weights for robust regression
#'
#' @param Y_data n x V data matrix
#' @param threshold Number of MADs for outlier detection
#' @param min_weight Minimum weight for outliers
#' @return n x V matrix of weights
#' @export
detect_outlier_timepoints <- function(Y_data, threshold = 3, min_weight = 0.1) {
  
  n <- nrow(Y_data)
  V <- ncol(Y_data)
  
  # Initialize weights
  weights <- matrix(1, n, V)
  
  # Detect outliers per voxel
  n_outliers_total <- 0
  
  for (v in 1:V) {
    y <- Y_data[, v]
    
    # Robust center and scale
    y_median <- median(y)
    y_mad <- median(abs(y - y_median)) * 1.4826
    
    if (y_mad > .Machine$double.eps) {
      # Standardized values
      z_scores <- abs(y - y_median) / y_mad
      
      # Identify outliers
      outliers <- which(z_scores > threshold)
      n_outliers_total <- n_outliers_total + length(outliers)
      
      # Soft downweighting using Huber-like weights
      if (length(outliers) > 0) {
        # Linear downweighting beyond threshold
        weights[outliers, v] <- pmax(min_weight, 
                                    1 - (z_scores[outliers] - threshold) / threshold)
      }
    }
  }
  
  if (n_outliers_total > 0) {
    message(sprintf("Detected %d outlier timepoints (%.1f%% of data)",
                   n_outliers_total, 100 * n_outliers_total / (n * V)))
  }
  
  return(weights)
}


#' Screen Voxels for Quality
#'
#' Identifies voxels that should be excluded or flagged
#'
#' @param Y_data n x V data matrix
#' @param min_variance Minimum temporal variance threshold
#' @param max_spike_fraction Maximum fraction of spike-like values
#' @return List with keep/flag indicators and quality metrics
#' @export
screen_voxels <- function(Y_data, 
                         min_variance = 1e-6,
                         max_spike_fraction = 0.1) {
  
  n <- nrow(Y_data)
  V <- ncol(Y_data)
  
  # Initialize
  keep_voxel <- rep(TRUE, V)
  flag_voxel <- rep(FALSE, V)
  quality_scores <- numeric(V)
  
  for (v in 1:V) {
    y <- Y_data[, v]
    
    # Handle NA/Inf values
    if (any(!is.finite(y))) {
      keep_voxel[v] <- FALSE
      quality_scores[v] <- 0
      flag_voxel[v] <- TRUE
      next
    }
    
    # Check variance
    y_var <- var(y)
    if (is.na(y_var) || y_var < min_variance) {
      keep_voxel[v] <- FALSE
      quality_scores[v] <- 0
      next
    }
    
    # Check for spikes (sudden large changes)
    y_diff <- abs(diff(y))
    y_diff_median <- median(y_diff)
    y_diff_mad <- median(abs(y_diff - y_diff_median)) * 1.4826
    
    if (y_diff_mad > 0) {
      spike_threshold <- y_diff_median + 5 * y_diff_mad
      n_spikes <- sum(y_diff > spike_threshold)
      spike_fraction <- n_spikes / (n - 1)
      
      if (spike_fraction > max_spike_fraction) {
        flag_voxel[v] <- TRUE
      }
    }
    
    # Quality score based on temporal smoothness
    quality_scores[v] <- 1 / (1 + mean(y_diff^2) / y_var)
  }
  
  n_excluded <- sum(!keep_voxel)
  n_flagged <- sum(flag_voxel)
  
  if (n_excluded > 0) {
    warning(sprintf("Excluded %d low-quality voxels", n_excluded))
  }
  if (n_flagged > 0) {
    warning(sprintf("Flagged %d voxels with potential artifacts", n_flagged))
  }
  
  return(list(
    keep = keep_voxel,
    flag = flag_voxel,
    quality_scores = quality_scores,
    n_excluded = n_excluded,
    n_flagged = n_flagged
  ))
}


#' Robust Regression with Huber Weights
#'
#' Implements robust regression using iteratively reweighted least squares
#'
#' @param X Design matrix
#' @param y Response vector
#' @param lambda Ridge parameter
#' @param huber_k Huber's k parameter
#' @param max_iter Maximum iterations
#' @return Regression coefficients
#' @keywords internal
robust_ridge_regression <- function(X, y, lambda = 0, huber_k = 1.345, max_iter = 10) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Initialize with standard ridge
  XtX <- crossprod(X)
  XtX_reg <- XtX + lambda * diag(p)
  beta <- solve(XtX_reg, crossprod(X, y))
  
  # Iteratively reweighted least squares
  for (iter in 1:max_iter) {
    # Compute residuals
    residuals <- y - X %*% beta
    
    # Robust scale estimate
    res_mad <- median(abs(residuals)) * 1.4826
    
    if (res_mad < .Machine$double.eps) break
    
    # Huber weights
    standardized_res <- abs(residuals) / res_mad
    weights <- ifelse(standardized_res <= huber_k, 
                     1, 
                     huber_k / standardized_res)
    
    # Weighted regression
    W <- diag(as.vector(weights))
    XtWX <- crossprod(X, W %*% X)
    XtWX_reg <- XtWX + lambda * diag(p)
    
    beta_new <- solve(XtWX_reg, crossprod(X, W %*% y))
    
    # Check convergence
    if (max(abs(beta_new - beta)) < 1e-6) break
    
    beta <- beta_new
  }
  
  return(as.vector(beta))
}
