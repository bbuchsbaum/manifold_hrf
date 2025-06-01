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
    cor_matrix <- matrix(0, N, N)
    max_cor <- 0
    n_near_duplicates <- 0
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
    if (!is.null(quality$has_duplicates) && quality$has_duplicates) {
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
#' @param preset Character string: "conservative", "balanced", "aggressive", 
#'   "fast", "quality", or "robust"
#' @param data_scale Optional data scale for parameter adaptation
#' @param n_voxels Optional number of voxels for memory-aware settings
#' @return List of parameters with metadata
#' @export
get_preset_params <- function(preset = c("conservative", "balanced", "aggressive",
                                        "fast", "quality", "robust"),
                             data_scale = NULL,
                             n_voxels = NULL) {
  
  preset <- match.arg(preset)
  
  base_params <- switch(preset,
    conservative = list(
      # Core parameters
      m_manifold_dim_target = 3,
      m_manifold_dim_min_variance = 0.90,
      lambda_gamma = 0.1,
      lambda_spatial_smooth = 1.0,
      lambda_beta_final = 0.01,
      lambda_ridge_Alss = 1e-5,
      k_local_nn_for_sigma = 10,
      num_neighbors_Lsp = 6,
      # Robustness settings
      use_robust_svd = TRUE,
      screen_voxels = TRUE,
      apply_hrf_constraints = TRUE,
      outlier_threshold = 3,
      max_iterations = 1,
      # Metadata
      description = "Conservative: Maximum stability, may underfit",
      use_case = "Noisy data, clinical studies, first-pass analysis"
    ),
    balanced = list(
      # Core parameters
      m_manifold_dim_target = 5,
      m_manifold_dim_min_variance = 0.95,
      lambda_gamma = 0.01,
      lambda_spatial_smooth = 0.5,
      lambda_beta_final = 0.001,
      lambda_ridge_Alss = 1e-6,
      k_local_nn_for_sigma = 7,
      num_neighbors_Lsp = 12,
      # Robustness settings
      use_robust_svd = TRUE,
      screen_voxels = TRUE,
      apply_hrf_constraints = TRUE,
      outlier_threshold = 3.5,
      max_iterations = 2,
      # Metadata
      description = "Balanced: Good tradeoff between stability and flexibility",
      use_case = "Standard fMRI studies, moderate noise levels"
    ),
    aggressive = list(
      # Core parameters
      m_manifold_dim_target = 8,
      m_manifold_dim_min_variance = 0.99,
      lambda_gamma = 0.001,
      lambda_spatial_smooth = 0.1,
      lambda_beta_final = 0.0001,
      lambda_ridge_Alss = 1e-7,
      k_local_nn_for_sigma = 5,
      num_neighbors_Lsp = 18,
      # Robustness settings
      use_robust_svd = FALSE,
      screen_voxels = TRUE,
      apply_hrf_constraints = FALSE,
      outlier_threshold = 4,
      max_iterations = 5,
      # Metadata
      description = "Aggressive: Maximum flexibility, requires clean data",
      use_case = "High-quality data, experienced users, fine-tuning"
    ),
    fast = list(
      # Core parameters
      m_manifold_dim_target = 4,
      m_manifold_dim_min_variance = 0.90,
      lambda_gamma = 0.05,
      lambda_spatial_smooth = 0.3,
      lambda_beta_final = 0.005,
      lambda_ridge_Alss = 1e-6,
      k_local_nn_for_sigma = 5,
      num_neighbors_Lsp = 6,
      # Robustness settings
      use_robust_svd = FALSE,
      screen_voxels = TRUE,
      apply_hrf_constraints = FALSE,
      outlier_threshold = 4,
      max_iterations = 1,
      # Performance settings
      use_parallel = TRUE,
      chunk_size = 5000,
      show_progress = TRUE,
      # Metadata
      description = "Fast: Optimized for speed, basic quality checks",
      use_case = "Large datasets, real-time analysis, exploratory work"
    ),
    quality = list(
      # Core parameters
      m_manifold_dim_target = 6,
      m_manifold_dim_min_variance = 0.97,
      lambda_gamma = 0.02,
      lambda_spatial_smooth = 0.7,
      lambda_beta_final = 0.002,
      lambda_ridge_Alss = 1e-6,
      k_local_nn_for_sigma = 10,
      num_neighbors_Lsp = 18,
      # Robustness settings
      use_robust_svd = TRUE,
      screen_voxels = TRUE,
      apply_hrf_constraints = TRUE,
      outlier_threshold = 3,
      max_iterations = 3,
      # Quality settings
      convergence_tol = 1e-5,
      min_voxel_var = 1e-6,
      hrf_peak_range = c(3, 9),
      # Metadata
      description = "Quality: Maximum accuracy, comprehensive checks",
      use_case = "Publication-quality results, small ROIs, detailed analysis"
    ),
    robust = list(
      # Core parameters
      m_manifold_dim_target = 4,
      m_manifold_dim_min_variance = 0.93,
      lambda_gamma = 0.05,
      lambda_spatial_smooth = 0.8,
      lambda_beta_final = 0.005,
      lambda_ridge_Alss = 1e-5,
      k_local_nn_for_sigma = 8,
      num_neighbors_Lsp = 10,
      # Robustness settings
      use_robust_svd = TRUE,
      screen_voxels = TRUE,
      apply_hrf_constraints = TRUE,
      outlier_threshold = 2.5,
      max_iterations = 2,
      # Fallback settings
      fallback_to_pca = TRUE,
      handle_zero_voxels = "noise",
      adaptive_smoothing = TRUE,
      # Metadata
      description = "Robust: Maximum reliability, handles difficult data",
      use_case = "Clinical data, motion artifacts, heterogeneous quality"
    )
  )
  
  # Scale parameters if data scale provided
  if (!is.null(data_scale) && data_scale > 0) {
    scale_factor <- data_scale^2
    base_params$lambda_gamma <- base_params$lambda_gamma * scale_factor
    base_params$lambda_beta_final <- base_params$lambda_beta_final * scale_factor
    base_params$lambda_ridge_Alss <- base_params$lambda_ridge_Alss * scale_factor
  }
  
  # Adjust for data size if provided
  if (!is.null(n_voxels) && n_voxels > 0) {
    if (n_voxels > 50000) {
      base_params$chunk_size <- 2000
      base_params$use_parallel <- TRUE
      message("Large dataset detected: Enabling chunked processing")
    }
    if (n_voxels < 1000) {
      base_params$num_neighbors_Lsp <- min(base_params$num_neighbors_Lsp, 
                                           round(n_voxels / 10))
    }
  }
  
  # Add workflow functions
  base_params$print_summary <- function() {
    cat("\n=== M-HRF-LSS Parameter Preset ===\n")
    cat(sprintf("Preset: %s\n", preset))
    cat(sprintf("Description: %s\n", base_params$description))
    cat(sprintf("Use case: %s\n", base_params$use_case))
    cat("\nKey parameters:\n")
    cat(sprintf("  Manifold dimensions: %d (%.0f%% variance)\n", 
                base_params$m_manifold_dim_target,
                base_params$m_manifold_dim_min_variance * 100))
    cat(sprintf("  Regularization (gamma): %.3f\n", base_params$lambda_gamma))
    cat(sprintf("  Spatial smoothing: %.2f\n", base_params$lambda_spatial_smooth))
    cat(sprintf("  Robustness features: %s\n",
                ifelse(base_params$use_robust_svd, "Enabled", "Disabled")))
    cat("==================================\n\n")
  }
  
  base_params$validate_data <- function(Y_data, X_design_list) {
    cat("Validating data compatibility with preset...\n")
    
    n <- nrow(Y_data)
    V <- ncol(Y_data)
    
    issues <- character()
    
    # Check data size
    if (V > 50000 && preset == "quality") {
      issues <- c(issues, "Large dataset with 'quality' preset may be slow")
    }
    
    # Check for problematic voxels
    zero_check <- handle_zero_voxels(Y_data)
    if (zero_check$percent_problematic > 20) {
      issues <- c(issues, sprintf("%.1f%% voxels are zero/low-variance", 
                                 zero_check$percent_problematic))
    }
    
    # Check design rank
    for (i in seq_along(X_design_list)) {
      rank_check <- check_design_rank(X_design_list[[i]])
      if (rank_check$is_rank_deficient) {
        issues <- c(issues, sprintf("Design matrix %d is rank deficient", i))
      }
    }
    
    if (length(issues) > 0) {
      cat("⚠️  Potential issues detected:\n")
      for (issue in issues) {
        cat(sprintf("   - %s\n", issue))
      }
      cat("\nConsider using 'robust' preset or addressing these issues.\n")
    } else {
      cat("✓ Data appears compatible with preset\n")
    }
    
    invisible(list(compatible = length(issues) == 0, issues = issues))
  }
  
  message(sprintf("Using '%s' parameter preset", preset))
  if (!is.null(base_params$description)) {
    message(base_params$description)
  }
  
  class(base_params) <- c("mhrf_preset", "list")
  return(base_params)
}


#' Print method for mhrf_preset
#' 
#' @param x mhrf_preset object
#' @param ... Additional arguments (ignored)
#' @export
print.mhrf_preset <- function(x, ...) {
  if (is.function(x$print_summary)) {
    x$print_summary()
  } else {
    # Fallback printing
    cat("M-HRF-LSS Parameter Preset\n")
    str(x[!sapply(x, is.function)])
  }
  invisible(x)
}


#' Create Full Workflow Configuration
#'
#' Creates a complete workflow configuration with all settings
#'
#' @param preset Base preset to use
#' @param custom_params List of custom parameter overrides
#' @param data_checks Whether to enable data validation checks
#' @param output_dir Directory for outputs (QC reports, logs, etc.)
#' @return Complete workflow configuration
#' @export
create_workflow_config <- function(preset = "balanced",
                                 custom_params = NULL,
                                 data_checks = TRUE,
                                 output_dir = NULL) {
  
  # Start with preset
  config <- get_preset_params(preset)
  
  # Override with custom parameters
  if (!is.null(custom_params)) {
    for (param in names(custom_params)) {
      config[[param]] <- custom_params[[param]]
    }
  }
  
  # Add workflow settings
  config$data_checks <- data_checks
  config$output_dir <- output_dir %||% file.path(getwd(), "mhrf_output")
  config$save_intermediates <- FALSE
  config$generate_qc_report <- TRUE
  config$timestamp <- Sys.time()
  
  # Create output directory if needed
  if (!dir.exists(config$output_dir)) {
    dir.create(config$output_dir, recursive = TRUE)
    message(sprintf("Created output directory: %s", config$output_dir))
  }
  
  # Add logging function
  config$log_file <- file.path(config$output_dir, 
                              sprintf("mhrf_log_%s.txt", 
                                     format(config$timestamp, "%Y%m%d_%H%M%S")))
  
  config$log_message <- function(msg, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    log_entry <- sprintf("[%s] %s: %s\n", timestamp, level, msg)
    cat(log_entry)
    cat(log_entry, file = config$log_file, append = TRUE)
  }
  
  # Initialize log
  config$log_message(sprintf("M-HRF-LSS workflow initialized with '%s' preset", preset))
  
  class(config) <- c("mhrf_workflow_config", class(config))
  return(config)
}


#' Apply Physiological HRF Constraints
#'
#' Ensures HRF shapes are physiologically plausible
#'
#' @param hrf_matrix p x V matrix of HRF shapes (timepoints x voxels)
#' @param TR Time repetition in seconds
#' @param peak_range Acceptable peak time range in seconds (default c(2, 10))
#' @param enforce_positive Whether to enforce positive post-peak values
#' @param project_to_plausible Whether to project to plausible subspace
#' @return List with constrained HRFs and quality metrics
#' @export
apply_hrf_physiological_constraints <- function(hrf_matrix,
                                               TR = 2,
                                               peak_range = c(2, 10),
                                               enforce_positive = TRUE,
                                               project_to_plausible = TRUE) {
  
  p <- nrow(hrf_matrix)
  V <- ncol(hrf_matrix)
  time_vec <- (0:(p-1)) * TR
  
  # Initialize output
  hrf_constrained <- hrf_matrix
  quality_metrics <- matrix(NA, 4, V)
  rownames(quality_metrics) <- c("peak_time", "is_plausible", "adjustment_made", "integral")
  
  for (v in 1:V) {
    hrf <- hrf_matrix[, v]
    
    # Skip zero HRFs
    if (all(abs(hrf) < .Machine$double.eps)) {
      quality_metrics["is_plausible", v] <- FALSE
      next
    }
    
    # Find peak
    peak_idx <- which.max(hrf)
    peak_time <- time_vec[peak_idx]
    quality_metrics["peak_time", v] <- peak_time
    
    # Check if peak is in reasonable range
    peak_ok <- peak_time >= peak_range[1] && peak_time <= peak_range[2]
    
    # Check if integral is positive (net positive response)
    integral <- sum(hrf)
    quality_metrics["integral", v] <- integral
    integral_ok <- integral > 0
    
    # Check post-peak positivity (after initial undershoot)
    post_peak_ok <- TRUE
    if (peak_idx < p - 5) {  # Need at least 5 points after peak
      # Allow for initial undershoot but check late response
      late_response <- hrf[(peak_idx + 5):p]
      post_peak_ok <- mean(late_response) > -0.1 * max(hrf)
    }
    
    quality_metrics["is_plausible", v] <- peak_ok && integral_ok && post_peak_ok
    quality_metrics["adjustment_made", v] <- 0
    
    # Apply corrections if needed
    if (project_to_plausible && !quality_metrics["is_plausible", v]) {
      
      # Correct peak time if needed
      if (!peak_ok) {
        if (peak_time < peak_range[1]) {
          # Shift HRF to right
          shift_samples <- ceiling((peak_range[1] - peak_time) / TR)
          hrf_new <- c(rep(0, shift_samples), hrf[1:(p - shift_samples)])
        } else {
          # Shift HRF to left
          shift_samples <- ceiling((peak_time - peak_range[2]) / TR)
          hrf_new <- c(hrf[(shift_samples + 1):p], rep(0, shift_samples))
        }
        hrf_constrained[, v] <- hrf_new
        quality_metrics["adjustment_made", v] <- 1
      }
      
      # Ensure positive integral
      if (!integral_ok) {
        # Add small positive offset
        hrf_constrained[, v] <- hrf_constrained[, v] + 0.01
        quality_metrics["adjustment_made", v] <- 1
      }
      
      # Fix negative late response
      if (enforce_positive && !post_peak_ok && peak_idx < p - 5) {
        late_idx <- (peak_idx + 5):p
        hrf_constrained[late_idx, v] <- pmax(hrf_constrained[late_idx, v], 0)
        quality_metrics["adjustment_made", v] <- 1
      }
    }
  }
  
  # Compute overall reasonableness score
  n_plausible <- sum(quality_metrics["is_plausible", ], na.rm = TRUE)
  n_adjusted <- sum(quality_metrics["adjustment_made", ], na.rm = TRUE)
  
  return(list(
    hrf_constrained = hrf_constrained,
    quality_metrics = quality_metrics,
    percent_plausible = 100 * n_plausible / V,
    percent_adjusted = 100 * n_adjusted / V
  ))
}


#' Compute HRF Reasonableness Score
#'
#' Computes a score indicating how physiologically reasonable an HRF is
#'
#' @param hrf_vector Single HRF time course
#' @param TR Time repetition in seconds
#' @return Reasonableness score between 0 and 1
#' @keywords internal
compute_hrf_reasonableness <- function(hrf_vector, TR = 2) {
  
  p <- length(hrf_vector)
  time_vec <- (0:(p-1)) * TR
  
  # Zero HRF gets score 0
  if (all(abs(hrf_vector) < .Machine$double.eps)) {
    return(0)
  }
  
  # Normalize for comparison
  hrf_norm <- hrf_vector / max(abs(hrf_vector))
  
  # Score components
  scores <- numeric(5)
  
  # 1. Peak time score (best around 5s)
  peak_idx <- which.max(hrf_norm)
  peak_time <- time_vec[peak_idx]
  scores[1] <- exp(-0.5 * ((peak_time - 5) / 2)^2)  # Gaussian centered at 5s
  
  # 2. Peak width score (not too narrow or wide)
  half_max <- max(hrf_norm) / 2
  above_half <- which(hrf_norm > half_max)
  if (length(above_half) > 1) {
    fwhm <- (max(above_half) - min(above_half)) * TR
    scores[2] <- exp(-0.5 * ((fwhm - 6) / 3)^2)  # Ideal FWHM around 6s
  } else {
    scores[2] <- 0
  }
  
  # 3. Undershoot score (should have mild undershoot)
  if (peak_idx < p - 5) {
    undershoot_region <- hrf_norm[(peak_idx + 3):min(peak_idx + 10, p)]
    undershoot_depth <- -min(undershoot_region)
    # Ideal undershoot around 20% of peak
    scores[3] <- exp(-2 * abs(undershoot_depth - 0.2))
  } else {
    scores[3] <- 0.5  # Can't evaluate
  }
  
  # 4. Return to baseline score
  if (p > 15) {
    late_response <- mean(abs(hrf_norm[(p-5):p]))
    scores[4] <- exp(-10 * late_response)  # Should be near zero
  } else {
    scores[4] <- 0.5
  }
  
  # 5. Smoothness score (not too jagged)
  diff2 <- diff(diff(hrf_norm))
  roughness <- sum(diff2^2) / (p - 2)
  scores[5] <- exp(-5 * roughness)
  
  # Weighted average
  weights <- c(0.3, 0.2, 0.2, 0.15, 0.15)
  overall_score <- sum(scores * weights)
  
  return(overall_score)
}


#' Track Convergence Metrics
#'
#' Tracks convergence metrics across iterations for diagnostics
#'
#' @param current_values Current parameter values (vector or matrix)
#' @param previous_values Previous parameter values
#' @param iteration Current iteration number
#' @param metric_name Name of the metric being tracked
#' @param history Optional existing history to append to
#' @return Updated convergence history
#' @export
track_convergence_metrics <- function(current_values, 
                                    previous_values = NULL,
                                    iteration = 1,
                                    metric_name = "parameters",
                                    history = NULL) {
  
  # Initialize history if needed
  if (is.null(history)) {
    history <- list(
      iterations = numeric(),
      relative_change = numeric(),
      absolute_change = numeric(),
      max_change = numeric(),
      metric_name = metric_name,
      converged = FALSE,
      converged_at = NA
    )
  }
  
  history$iterations <- c(history$iterations, iteration)
  
  # Compute changes if previous values provided
  if (!is.null(previous_values)) {
    diff_vals <- current_values - previous_values
    
    # Relative change (avoid division by zero)
    denom <- pmax(abs(previous_values), 1e-10)
    rel_change <- sqrt(mean((diff_vals / denom)^2))
    
    # Absolute change
    abs_change <- sqrt(mean(diff_vals^2))
    
    # Maximum change
    max_change <- max(abs(diff_vals))
    
    history$relative_change <- c(history$relative_change, rel_change)
    history$absolute_change <- c(history$absolute_change, abs_change)
    history$max_change <- c(history$max_change, max_change)
  } else {
    # First iteration
    history$relative_change <- c(history$relative_change, NA)
    history$absolute_change <- c(history$absolute_change, NA)
    history$max_change <- c(history$max_change, NA)
  }
  
  return(history)
}


#' Check Convergence Status
#'
#' Checks if convergence criteria are met and provides diagnostics
#'
#' @param history Convergence history from track_convergence_metrics
#' @param rel_tol Relative tolerance for convergence
#' @param abs_tol Absolute tolerance for convergence
#' @param min_iterations Minimum iterations before checking convergence
#' @param patience Number of iterations to wait for improvement
#' @return List with convergence status and diagnostics
#' @export
check_convergence_status <- function(history,
                                   rel_tol = 1e-4,
                                   abs_tol = 1e-6,
                                   min_iterations = 2,
                                   patience = 3) {
  
  n_iter <- length(history$iterations)
  
  # Need minimum iterations
  if (n_iter < min_iterations) {
    return(list(
      converged = FALSE,
      reason = "insufficient_iterations",
      diagnostic = sprintf("Only %d iterations completed (min: %d)", 
                          n_iter, min_iterations)
    ))
  }
  
  # Get recent changes
  recent_rel <- tail(history$relative_change[!is.na(history$relative_change)], patience)
  recent_abs <- tail(history$absolute_change[!is.na(history$absolute_change)], patience)
  
  # Check relative tolerance
  if (length(recent_rel) >= patience && all(recent_rel < rel_tol)) {
    return(list(
      converged = TRUE,
      reason = "relative_tolerance",
      diagnostic = sprintf("Relative change < %.2e for %d iterations", 
                          rel_tol, patience),
      final_change = tail(recent_rel, 1)
    ))
  }
  
  # Check absolute tolerance
  if (length(recent_abs) >= patience && all(recent_abs < abs_tol)) {
    return(list(
      converged = TRUE,
      reason = "absolute_tolerance",
      diagnostic = sprintf("Absolute change < %.2e for %d iterations", 
                          abs_tol, patience),
      final_change = tail(recent_abs, 1)
    ))
  }
  
  # Check for stagnation (suspicious lack of change)
  if (length(recent_rel) >= 2) {
    rel_var <- var(recent_rel)
    if (rel_var < 1e-10) {
      return(list(
        converged = FALSE,
        reason = "stagnation",
        diagnostic = "Changes are suspiciously constant - may be stuck",
        variance = rel_var
      ))
    }
  }
  
  # Not converged yet
  return(list(
    converged = FALSE,
    reason = "in_progress",
    diagnostic = sprintf("Iteration %d: rel_change = %.2e, abs_change = %.2e",
                        n_iter, 
                        tail(recent_rel, 1),
                        tail(recent_abs, 1))
  ))
}


#' Compute Solution Quality Metrics
#'
#' Computes quality metrics for the current solution
#'
#' @param Y_data Data matrix (n x V)
#' @param Y_predicted Predicted data matrix (n x V)
#' @param Xi_matrix Manifold coordinates (m x V)
#' @param Beta_matrix Amplitude matrix (k x V)
#' @param lambda_smooth Smoothing parameter used
#' @return List of quality metrics
#' @export
compute_solution_quality <- function(Y_data, 
                                   Y_predicted,
                                   Xi_matrix = NULL,
                                   Beta_matrix = NULL,
                                   lambda_smooth = NULL) {
  
  metrics <- list()
  
  # Reconstruction error
  residuals <- Y_data - Y_predicted
  metrics$rmse <- sqrt(mean(residuals^2))
  metrics$r_squared <- 1 - sum(residuals^2) / sum((Y_data - mean(Y_data))^2)
  
  # Per-voxel R-squared
  V <- ncol(Y_data)
  metrics$r_squared_voxels <- numeric(V)
  for (v in 1:V) {
    ss_res <- sum(residuals[, v]^2)
    ss_tot <- sum((Y_data[, v] - mean(Y_data[, v]))^2)
    metrics$r_squared_voxels[v] <- 1 - ss_res / ss_tot
  }
  
  # Smoothness metrics if Xi provided
  if (!is.null(Xi_matrix)) {
    # Spatial smoothness (lower is smoother)
    m <- nrow(Xi_matrix)
    xi_smoothness <- numeric(m)
    for (i in 1:m) {
      # Simple measure: variance of spatial gradients
      xi_row <- Xi_matrix[i, ]
      if (length(xi_row) > 1) {
        xi_smoothness[i] <- var(diff(xi_row))
      }
    }
    metrics$xi_smoothness <- mean(xi_smoothness)
    
    # Effective degrees of freedom (complexity)
    # Approximated by ratio of variance explained
    metrics$effective_df <- sum(apply(Xi_matrix, 1, var) > 1e-6)
  }
  
  # Amplitude statistics if Beta provided
  if (!is.null(Beta_matrix)) {
    metrics$beta_mean <- mean(Beta_matrix)
    metrics$beta_sd <- sd(Beta_matrix)
    metrics$beta_sparsity <- mean(abs(Beta_matrix) < 0.1)
  }
  
  # Overall quality score (0-1)
  quality_components <- c(
    min(metrics$r_squared, 1),  # Fit quality
    1 - min(metrics$xi_smoothness / 10, 1),  # Smoothness (normalized)
    1 - metrics$beta_sparsity  # Non-sparsity
  )
  metrics$overall_quality <- mean(quality_components, na.rm = TRUE)
  
  return(metrics)
}


#' Check Design Matrix Rank
#'
#' Checks rank deficiency in design matrices and removes collinear columns
#'
#' @param X_design Design matrix to check
#' @param tol Tolerance for rank determination
#' @param remove_collinear Whether to remove collinear columns
#' @return List with checked/cleaned design matrix and diagnostics
#' @export
check_design_rank <- function(X_design, tol = 1e-10, remove_collinear = TRUE) {
  
  n <- nrow(X_design)
  p <- ncol(X_design)
  
  # Compute SVD for rank check
  svd_X <- svd(X_design)
  rank_X <- sum(svd_X$d > tol * svd_X$d[1])
  
  result <- list(
    original_rank = rank_X,
    full_rank = p,
    is_rank_deficient = rank_X < p,
    condition_number = svd_X$d[1] / svd_X$d[rank_X]
  )
  
  if (result$is_rank_deficient) {
    warning(sprintf("Design matrix is rank deficient: rank %d < %d columns", 
                   rank_X, p))
    
    if (remove_collinear) {
      # Identify which columns to keep
      # Use QR decomposition with pivoting
      qr_X <- qr(X_design, tol = tol)
      keep_cols <- qr_X$pivot[1:rank_X]
      
      result$X_cleaned <- X_design[, keep_cols, drop = FALSE]
      result$removed_columns <- setdiff(1:p, keep_cols)
      result$kept_columns <- keep_cols
      
      message(sprintf("Removed %d collinear columns: %s", 
                     length(result$removed_columns),
                     paste(result$removed_columns, collapse = ", ")))
    }
  } else {
    result$X_cleaned <- X_design
    result$removed_columns <- integer(0)
    result$kept_columns <- 1:p
  }
  
  return(result)
}


#' Handle Zero Voxels
#'
#' Identifies and handles zero-variance and all-zero voxels
#'
#' @param Y_data Data matrix (n x V)
#' @param min_variance Minimum variance threshold
#' @param replace_with How to handle zero voxels ("skip", "noise", or "mean")
#' @return List with cleaned data and zero voxel indices
#' @export
handle_zero_voxels <- function(Y_data, min_variance = 1e-8, replace_with = "skip") {
  
  n <- nrow(Y_data)
  V <- ncol(Y_data)
  
  # Identify problematic voxels
  voxel_vars <- apply(Y_data, 2, var)
  voxel_means <- colMeans(Y_data)
  
  all_zero <- abs(voxel_means) < .Machine$double.eps & voxel_vars < .Machine$double.eps
  low_variance <- voxel_vars < min_variance
  
  zero_indices <- which(all_zero | low_variance)
  n_zero <- sum(all_zero)
  n_low_var <- sum(low_variance & !all_zero)
  
  # Report findings
  if (length(zero_indices) > 0) {
    message(sprintf("Found %d problematic voxels: %d all-zero, %d low-variance", 
                   length(zero_indices), n_zero, n_low_var))
  }
  
  # Handle based on strategy
  Y_cleaned <- Y_data
  if (length(zero_indices) > 0 && replace_with != "skip") {
    if (replace_with == "noise") {
      # Replace with small noise
      good_vars <- voxel_vars[voxel_vars > min_variance & is.finite(voxel_vars)]
      if (length(good_vars) > 0) {
        noise_sd <- sqrt(median(good_vars) * 0.01)
      } else {
        noise_sd <- 0.1  # Fallback value
      }
      Y_cleaned[, zero_indices] <- matrix(rnorm(n * length(zero_indices), sd = noise_sd),
                                         n, length(zero_indices))
    } else if (replace_with == "mean") {
      # Replace with global mean signal
      if (length(zero_indices) < V) {
        mean_signal <- rowMeans(Y_data[, -zero_indices, drop = FALSE], na.rm = TRUE)
        Y_cleaned[, zero_indices] <- mean_signal
      }
    }
  }
  
  return(list(
    Y_cleaned = Y_cleaned,
    zero_indices = zero_indices,
    all_zero_indices = which(all_zero),
    low_var_indices = which(low_variance & !all_zero),
    n_problematic = length(zero_indices),
    percent_problematic = 100 * length(zero_indices) / V
  ))
}


#' Monitor Condition Numbers
#'
#' Monitors condition numbers throughout the pipeline
#'
#' @param matrix_list Named list of matrices to check
#' @param warn_threshold Threshold for warning (default 1e8)
#' @param error_threshold Threshold for error (default 1e12)
#' @return Data frame with condition monitoring results
#' @export
monitor_condition_numbers <- function(matrix_list, 
                                    warn_threshold = 1e8,
                                    error_threshold = 1e12) {
  
  results <- data.frame(
    matrix_name = character(),
    dimensions = character(),
    condition_number = numeric(),
    status = character(),
    stringsAsFactors = FALSE
  )
  
  for (name in names(matrix_list)) {
    mat <- matrix_list[[name]]
    
    if (is.matrix(mat) || inherits(mat, "Matrix")) {
      # Compute condition number
      kappa_val <- tryCatch({
        if (nrow(mat) == ncol(mat)) {
          # Square matrix - use standard condition number
          kappa(mat)
        } else {
          # Rectangular - use SVD-based condition
          sv <- svd(mat, nu = 0, nv = 0)
          sv$d[1] / sv$d[length(sv$d)]
        }
      }, error = function(e) NA)
      
      # Determine status
      status <- "OK"
      if (is.na(kappa_val)) {
        status <- "ERROR"
      } else if (kappa_val > error_threshold) {
        status <- "CRITICAL"
        warning(sprintf("Matrix '%s' has critical condition number: %.2e", 
                       name, kappa_val))
      } else if (kappa_val > warn_threshold) {
        status <- "WARNING"
        message(sprintf("Matrix '%s' has high condition number: %.2e", 
                       name, kappa_val))
      }
      
      results <- rbind(results, data.frame(
        matrix_name = name,
        dimensions = sprintf("%dx%d", nrow(mat), ncol(mat)),
        condition_number = kappa_val,
        status = status,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Summary message
  n_warning <- sum(results$status == "WARNING")
  n_critical <- sum(results$status == "CRITICAL")
  
  if (n_critical > 0) {
    warning(sprintf("Found %d matrices with critical conditioning!", n_critical))
  }
  if (n_warning > 0) {
    message(sprintf("Found %d matrices with poor conditioning", n_warning))
  }
  
  return(results)
}


#' Create Progress Bar
#'
#' Creates a simple progress bar for long operations
#'
#' @param total Total number of items
#' @param width Width of progress bar
#' @return Progress bar object
#' @export
create_progress_bar <- function(total, width = 50) {
  
  pb <- list(
    total = total,
    current = 0,
    width = width,
    start_time = Sys.time(),
    last_update = Sys.time()
  )
  
  class(pb) <- "simple_progress_bar"
  return(pb)
}


#' Update Progress Bar
#'
#' Updates and displays progress bar
#'
#' @param pb Progress bar object
#' @param increment Increment amount (default 1)
#' @param message Optional message to display
#' @return Updated progress bar (invisibly)
#' @export
update_progress_bar <- function(pb, increment = 1, message = NULL) {
  
  pb$current <- pb$current + increment
  
  # Only update display every 0.1 seconds to avoid spam
  if (difftime(Sys.time(), pb$last_update, units = "secs") < 0.1 && 
      pb$current < pb$total) {
    return(invisible(pb))
  }
  
  pb$last_update <- Sys.time()
  
  # Calculate progress
  progress <- pb$current / pb$total
  filled <- round(progress * pb$width)
  
  # Time estimation
  elapsed <- difftime(Sys.time(), pb$start_time, units = "secs")
  if (pb$current > 0) {
    eta <- elapsed * (pb$total / pb$current - 1)
    eta_str <- sprintf("ETA: %s", format_time(eta))
  } else {
    eta_str <- "ETA: --:--"
  }
  
  # Build progress bar
  bar <- sprintf("[%s%s] %d/%d (%.1f%%) %s",
                paste(rep("=", filled), collapse = ""),
                paste(rep(" ", pb$width - filled), collapse = ""),
                pb$current, pb$total,
                progress * 100,
                eta_str)
  
  # Add message if provided
  if (!is.null(message)) {
    bar <- paste(bar, "-", message)
  }
  
  # Update display
  cat("\r", bar, sep = "")
  
  # New line when complete
  if (pb$current >= pb$total) {
    cat("\n")
    total_time <- difftime(Sys.time(), pb$start_time, units = "secs")
    message(sprintf("Completed in %s", format_time(total_time)))
  }
  
  return(invisible(pb))
}


#' Format Time Helper
#'
#' Formats seconds into human-readable time
#'
#' @param seconds Number of seconds
#' @return Formatted time string
#' @keywords internal
format_time <- function(seconds) {
  if (seconds < 60) {
    return(sprintf("%.0fs", seconds))
  } else if (seconds < 3600) {
    return(sprintf("%.0fm %.0fs", seconds %/% 60, seconds %% 60))
  } else {
    return(sprintf("%.0fh %.0fm", seconds %/% 3600, (seconds %% 3600) %/% 60))
  }
}