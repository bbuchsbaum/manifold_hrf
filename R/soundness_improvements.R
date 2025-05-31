# Soundness Improvements for M-HRF-LSS
# Implementation of SOUND-* improvements

#' Check HRF Library Quality
#'
#' Evaluates the quality and diversity of an HRF library matrix
#'
#' @param L_library_matrix p x N matrix of HRF shapes
#' @return List with quality metrics and flags
#' @keywords internal
check_hrf_library_quality <- function(L_library_matrix) {
  
  p <- nrow(L_library_matrix)
  N <- ncol(L_library_matrix)
  
  quality <- list()
  
  # Check for duplicate HRFs
  # Handle case where all HRFs are constant
  tryCatch({
    cor_matrix <- cor(L_library_matrix)
    diag(cor_matrix) <- 0
    max_cor <- max(abs(cor_matrix))
    n_near_duplicates <- sum(abs(cor_matrix) > 0.99) / 2
  }, error = function(e) {
    # If correlation fails (e.g., zero variance), treat as degenerate
    cor_matrix <<- matrix(0, N, N)
    max_cor <<- 0
    n_near_duplicates <<- 0
  })
  
  quality$has_duplicates <- n_near_duplicates > 0
  quality$n_duplicates <- n_near_duplicates
  quality$max_correlation <- max_cor
  
  # Check condition number
  svd_L <- svd(L_library_matrix)
  condition_number <- svd_L$d[1] / svd_L$d[min(p, N)]
  quality$condition_number <- condition_number
  quality$is_ill_conditioned <- condition_number > 1e10
  
  # Check for all-zero or constant HRFs
  zero_hrfs <- apply(L_library_matrix, 2, function(x) all(x == 0))
  constant_hrfs <- apply(L_library_matrix, 2, function(x) sd(x) < .Machine$double.eps)
  
  quality$n_zero_hrfs <- sum(zero_hrfs)
  quality$n_constant_hrfs <- sum(constant_hrfs)
  quality$has_degenerate_hrfs <- any(zero_hrfs | constant_hrfs)
  
  # Check diversity (mean pairwise distance)
  if (N > 1) {
    distances <- dist(t(L_library_matrix))
    quality$mean_diversity <- mean(distances)
    quality$min_diversity <- min(distances)
    quality$is_low_diversity <- quality$min_diversity < 0.01
  } else {
    quality$mean_diversity <- NA
    quality$min_diversity <- NA
    quality$is_low_diversity <- TRUE
  }
  
  # Overall quality flag
  quality$is_good_quality <- !quality$has_duplicates && 
                            !quality$is_ill_conditioned && 
                            !quality$has_degenerate_hrfs && 
                            !quality$is_low_diversity
  
  return(quality)
}


#' Remove Duplicate HRFs from Library
#'
#' Removes near-duplicate HRFs based on correlation threshold
#'
#' @param L_library_matrix p x N matrix of HRF shapes
#' @param cor_threshold Correlation threshold for duplicates (default 0.99)
#' @return Cleaned matrix with duplicates removed
#' @keywords internal
remove_duplicate_hrfs <- function(L_library_matrix, cor_threshold = 0.99) {
  
  N <- ncol(L_library_matrix)
  if (N <= 1) return(L_library_matrix)
  
  # Compute correlations
  cor_matrix <- cor(L_library_matrix)
  
  # Find which HRFs to keep
  keep <- rep(TRUE, N)
  
  for (i in 1:(N-1)) {
    if (keep[i]) {
      # Mark duplicates of HRF i for removal
      duplicates <- which(cor_matrix[i, (i+1):N] > cor_threshold) + i
      if (length(duplicates) > 0) {
        keep[duplicates] <- FALSE
      }
    }
  }
  
  n_removed <- sum(!keep)
  if (n_removed > 0) {
    message(sprintf("Removed %d duplicate HRFs (correlation > %.2f)", 
                   n_removed, cor_threshold))
  }
  
  return(L_library_matrix[, keep, drop = FALSE])
}


#' Compute PCA Basis as Fallback
#'
#' Computes PCA-based reconstructor when manifold construction fails
#'
#' @param L_library_matrix p x N matrix of HRF shapes
#' @param m_target Target dimensionality
#' @param min_variance Minimum variance to retain
#' @return List with B_reconstructor_matrix and metadata
#' @keywords internal
compute_pca_fallback <- function(L_library_matrix, m_target, min_variance = 0.95) {
  
  message("Using PCA fallback for manifold construction")
  
  p <- nrow(L_library_matrix)
  N <- ncol(L_library_matrix)
  
  # Center the HRFs
  L_centered <- L_library_matrix - rowMeans(L_library_matrix)
  
  # Compute SVD
  svd_result <- svd(L_centered)
  
  # Determine dimensionality based on variance
  var_explained <- svd_result$d^2 / sum(svd_result$d^2)
  cum_var <- cumsum(var_explained)
  m_auto <- which(cum_var >= min_variance)[1]
  
  # Use minimum of target and auto-selected
  m_final <- min(m_target, m_auto, length(svd_result$d))
  
  # Ensure we explain at least min_variance
  if (cum_var[m_final] < min_variance && m_final < length(svd_result$d)) {
    m_final <- m_auto
  }
  
  # Extract components
  B_reconstructor <- svd_result$u[, 1:m_final, drop = FALSE]
  
  # Compute coordinates (for compatibility with manifold output)
  Phi_coords <- svd_result$v[, 1:m_final, drop = FALSE] %*% 
                diag(svd_result$d[1:m_final], m_final, m_final)
  
  return(list(
    B_reconstructor_matrix = B_reconstructor,
    Phi_coords_matrix = Phi_coords,
    eigenvalues_S_vector = c(1, var_explained),  # Add trivial eigenvalue
    m_final_dim = m_final,
    m_auto_selected_dim = m_auto,
    method_used = "PCA",
    variance_explained = cum_var[m_final]
  ))
}


#' Enhanced Manifold Basis Reconstructor with Fallback
#'
#' Robust version of get_manifold_basis_reconstructor_core with PCA fallback
#'
#' @param S_markov_matrix N x N Markov transition matrix
#' @param L_library_matrix p x N HRF library matrix
#' @param m_manifold_dim_target Target manifold dimensionality
#' @param m_manifold_dim_min_variance Minimum variance threshold
#' @param fallback_to_pca Whether to fall back to PCA on failure
#' @return List with reconstructor matrix and metadata
#' @export
get_manifold_basis_reconstructor_robust <- function(S_markov_matrix,
                                                   L_library_matrix,
                                                   m_manifold_dim_target,
                                                   m_manifold_dim_min_variance = 0.95,
                                                   fallback_to_pca = TRUE) {
  
  # First check library quality
  quality <- check_hrf_library_quality(L_library_matrix)
  
  if (!quality$is_good_quality) {
    warning("HRF library has quality issues:")
    if (quality$has_duplicates) {
      warning(sprintf("  - Found %d near-duplicate HRFs", quality$n_duplicates))
    }
    if (quality$is_ill_conditioned) {
      warning(sprintf("  - Library is ill-conditioned (condition number: %.2e)", 
                     quality$condition_number))
    }
    if (quality$has_degenerate_hrfs) {
      warning(sprintf("  - Found %d zero and %d constant HRFs", 
                     quality$n_zero_hrfs, quality$n_constant_hrfs))
    }
    if (quality$is_low_diversity) {
      warning("  - Low diversity in HRF shapes")
    }
  }
  
  # Try standard manifold construction
  result <- tryCatch({
    
    # Check if S_markov is degenerate
    if (is.matrix(S_markov_matrix)) {
      S_condition <- kappa(S_markov_matrix)
    } else {
      # For sparse matrices, check a subset
      S_sample <- as.matrix(S_markov_matrix[1:min(10, nrow(S_markov_matrix)), 
                                            1:min(10, ncol(S_markov_matrix))])
      S_condition <- kappa(S_sample)
    }
    
    if (S_condition > 1e12) {
      warning(sprintf("Markov matrix is poorly conditioned (kappa = %.2e)", S_condition))
      if (fallback_to_pca) {
        stop("Triggering PCA fallback due to poor conditioning")
      }
    }
    
    # Call original function
    result <- get_manifold_basis_reconstructor_core(
      S_markov_matrix = S_markov_matrix,
      L_library_matrix = L_library_matrix,
      m_manifold_dim_target = m_manifold_dim_target,
      m_manifold_dim_min_variance = m_manifold_dim_min_variance
    )
    
    # Check if result is reasonable
    if (any(is.na(result$B_reconstructor_matrix)) || 
        any(is.infinite(result$B_reconstructor_matrix))) {
      stop("Manifold construction produced invalid results")
    }
    
    # Add quality metrics
    result$library_quality <- quality
    result$method_used <- "diffusion_map"
    
    result
    
  }, error = function(e) {
    if (fallback_to_pca) {
      warning(sprintf("Manifold construction failed: %s", e$message))
      warning("Falling back to PCA-based approach")
      
      # Use PCA fallback
      pca_result <- compute_pca_fallback(
        L_library_matrix = L_library_matrix,
        m_target = m_manifold_dim_target,
        min_variance = m_manifold_dim_min_variance
      )
      
      pca_result$library_quality <- quality
      pca_result$original_error <- e$message
      
      return(pca_result)
    } else {
      stop(e)
    }
  })
  
  return(result)
}


#' Adaptive Parameter Selection
#'
#' Suggests parameters based on data characteristics
#'
#' @param Y_data n x V data matrix
#' @param X_design_list List of design matrices
#' @param voxel_coords V x 3 coordinate matrix (optional)
#' @return List of suggested parameters
#' @export
suggest_parameters <- function(Y_data, X_design_list = NULL, voxel_coords = NULL) {
  
  n <- nrow(Y_data)
  V <- ncol(Y_data)
  
  params <- list()
  
  # Estimate noise level from data
  # Use median absolute deviation for robustness
  Y_centered <- Y_data - rowMeans(Y_data)
  data_mad <- median(abs(Y_centered))
  data_scale <- data_mad * 1.4826  # MAD to SD conversion
  
  # Lambda for gamma regression - scale with noise
  params$lambda_gamma <- 0.01 * data_scale^2
  
  # Spatial smoothing - adapt to voxel density
  if (!is.null(voxel_coords) && V > 10) {
    # Estimate typical voxel spacing
    distances <- as.matrix(dist(voxel_coords))
    diag(distances) <- NA
    median_nn_dist <- median(apply(distances, 1, min, na.rm = TRUE))
    
    # Less smoothing for denser sampling
    if (median_nn_dist < 2) {
      params$lambda_spatial_smooth <- 0.1
    } else if (median_nn_dist < 4) {
      params$lambda_spatial_smooth <- 0.5
    } else {
      params$lambda_spatial_smooth <- 1.0
    }
    
    # Adjust neighbors based on density
    params$num_neighbors_Lsp <- min(26, max(6, round(V^(1/3))))
  } else {
    params$lambda_spatial_smooth <- 0.5
    params$num_neighbors_Lsp <- 6
  }
  
  # Final beta regularization
  params$lambda_beta_final <- params$lambda_gamma * 0.1
  
  # Ridge for LSS
  params$lambda_ridge_Alss <- 1e-6 * data_scale^2
  
  # Manifold dimension - based on data complexity
  if (!is.null(X_design_list)) {
    n_conditions <- length(X_design_list)
    # More conditions might need more manifold dimensions
    params$m_manifold_dim_target <- min(8, max(3, n_conditions + 2))
  } else {
    params$m_manifold_dim_target <- 5
  }
  
  # Report suggestions
  message("Suggested parameters based on data characteristics:")
  message(sprintf("  Data scale (MAD): %.3f", data_scale))
  message(sprintf("  lambda_gamma: %.4f", params$lambda_gamma))
  message(sprintf("  lambda_spatial_smooth: %.2f", params$lambda_spatial_smooth))
  message(sprintf("  lambda_beta_final: %.4f", params$lambda_beta_final))
  message(sprintf("  m_manifold_dim_target: %d", params$m_manifold_dim_target))
  
  return(params)
}


#' Preset Parameter Configurations
#'
#' Returns preset parameter configurations for different analysis styles
#'
#' @param preset Character string: "conservative", "balanced", or "aggressive"
#' @param data_scale Optional data scale for parameter adaptation
#' @return List of parameters
#' @export
get_preset_params <- function(preset = c("conservative", "balanced", "aggressive"),
                             data_scale = NULL) {
  
  preset <- match.arg(preset)
  
  base_params <- switch(preset,
    conservative = list(
      m_manifold_dim_target = 3,
      m_manifold_dim_min_variance = 0.90,
      lambda_gamma = 0.1,
      lambda_spatial_smooth = 1.0,
      lambda_beta_final = 0.01,
      lambda_ridge_Alss = 1e-5,
      k_local_nn_for_sigma = 10,
      num_neighbors_Lsp = 6
    ),
    balanced = list(
      m_manifold_dim_target = 5,
      m_manifold_dim_min_variance = 0.95,
      lambda_gamma = 0.01,
      lambda_spatial_smooth = 0.5,
      lambda_beta_final = 0.001,
      lambda_ridge_Alss = 1e-6,
      k_local_nn_for_sigma = 7,
      num_neighbors_Lsp = 12
    ),
    aggressive = list(
      m_manifold_dim_target = 8,
      m_manifold_dim_min_variance = 0.99,
      lambda_gamma = 0.001,
      lambda_spatial_smooth = 0.1,
      lambda_beta_final = 0.0001,
      lambda_ridge_Alss = 1e-7,
      k_local_nn_for_sigma = 5,
      num_neighbors_Lsp = 18
    )
  )
  
  # Scale parameters if data scale provided
  if (!is.null(data_scale) && data_scale > 0) {
    scale_factor <- data_scale^2
    base_params$lambda_gamma <- base_params$lambda_gamma * scale_factor
    base_params$lambda_beta_final <- base_params$lambda_beta_final * scale_factor
    base_params$lambda_ridge_Alss <- base_params$lambda_ridge_Alss * scale_factor
  }
  
  message(sprintf("Using '%s' parameter preset", preset))
  
  return(base_params)
}