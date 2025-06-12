# Configuration Helper Functions

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

  # Default sign alignment method
  base_params$ident_sign_method <- "canonical_correlation"

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


