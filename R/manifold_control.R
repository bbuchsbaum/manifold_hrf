#' Control Parameters for Manifold HRF Estimation
#'
#' @description
#' Constructs a list of control parameters for the manifold HRF estimation algorithm.
#' Provides preset configurations and allows fine-grained control over all parameters.
#'
#' @param preset Character string specifying a preset configuration:
#'   \describe{
#'     \item{"fast"}{Quick estimation with reduced iterations and neighbors}
#'     \item{"balanced"}{Default settings balancing speed and accuracy}
#'     \item{"thorough"}{High-quality estimation with more iterations}
#'   }
#' @param m_manifold_dim Integer; dimension of the HRF manifold (NULL for automatic selection)
#' @param lambda_manifold Numeric; regularization for manifold construction
#' @param num_neighbors_mfd Integer; number of neighbors for manifold construction
#' @param lambda_spatial_smooth Numeric; spatial smoothing regularization strength
#' @param num_neighbors_Lsp Integer; number of neighbors for spatial graph
#' @param lambda_gamma Numeric; ridge penalty for voxel-wise fitting
#' @param lambda_condition_betas Numeric; regularization for condition-level betas
#' @param max_iter Integer; maximum iterations for optimization
#' @param tolerance Numeric; convergence tolerance
#' @param orthogonal_approx Logical; use orthogonal approximation in GLM
#' @param n_cores Integer; number of cores for parallel processing
#' @param verbose_level Integer; verbosity level (0 = silent, 1 = progress, 2 = detailed)
#' @param ram_threshold_GB Numeric; RAM threshold for memory strategy selection
#' @param generate_qc_report Logical; whether to generate quality control metrics
#' @param ... Additional parameters to override preset values
#'
#' @return A list of class "manifold_control" containing all control parameters
#'
#' @examples
#' # Use balanced preset
#' ctrl <- manifold_control()
#' 
#' # Fast estimation
#' ctrl_fast <- manifold_control(preset = "fast")
#' 
#' # Custom settings
#' ctrl_custom <- manifold_control(
#'   preset = "balanced",
#'   max_iter = 200,
#'   lambda_spatial_smooth = 1.0
#' )
#'
#' @export
manifold_control <- function(preset = c("balanced", "fast", "thorough"),
                           m_manifold_dim = NULL,
                           lambda_manifold = 0.1,
                           num_neighbors_mfd = 75,
                           lambda_spatial_smooth = 0.5,
                           num_neighbors_Lsp = 6,
                           lambda_gamma = 1e-6,
                           lambda_condition_betas = 1e-4,
                           max_iter = 100,
                           tolerance = 1e-6,
                           orthogonal_approx = FALSE,
                           n_cores = 1,
                           verbose_level = 1,
                           ram_threshold_GB = 8,
                           generate_qc_report = FALSE,
                           ...) {
  
  preset <- match.arg(preset)
  
  # Define preset configurations
  presets <- list(
    fast = list(
      max_iter = 50,
      num_neighbors_mfd = 50,
      num_neighbors_Lsp = 4,
      tolerance = 1e-4,
      lambda_manifold = 0.2,
      lambda_spatial_smooth = 0.3
    ),
    balanced = list(
      max_iter = 100,
      num_neighbors_mfd = 75,
      num_neighbors_Lsp = 6,
      tolerance = 1e-6,
      lambda_manifold = 0.1,
      lambda_spatial_smooth = 0.5
    ),
    thorough = list(
      max_iter = 500,
      num_neighbors_mfd = 100,
      num_neighbors_Lsp = 8,
      tolerance = 1e-8,
      lambda_manifold = 0.05,
      lambda_spatial_smooth = 0.7
    )
  )
  
  # Start with preset values
  params <- presets[[preset]]
  
  # Override with explicitly provided arguments
  call_args <- as.list(match.call())[-1]
  call_args$preset <- NULL  # Remove preset from override list
  
  # Handle ... arguments
  extra_args <- list(...)
  
  # Combine all arguments
  all_args <- c(call_args, extra_args)
  
  # Override preset values with user-specified values
  for (arg_name in names(all_args)) {
    if (arg_name != "...") {
      params[[arg_name]] <- all_args[[arg_name]]
    }
  }
  
  # Add arguments that weren't in preset
  default_args <- list(
    m_manifold_dim = m_manifold_dim,
    lambda_gamma = lambda_gamma,
    lambda_condition_betas = lambda_condition_betas,
    orthogonal_approx = orthogonal_approx,
    n_cores = n_cores,
    verbose_level = verbose_level,
    ram_threshold_GB = ram_threshold_GB,
    generate_qc_report = generate_qc_report
  )
  
  for (arg_name in names(default_args)) {
    if (!(arg_name %in% names(params))) {
      params[[arg_name]] <- default_args[[arg_name]]
    }
  }
  
  # Validation
  if (!is.null(params$m_manifold_dim)) {
    if (!is.numeric(params$m_manifold_dim) || params$m_manifold_dim < 1) {
      stop("m_manifold_dim must be a positive integer or NULL")
    }
    params$m_manifold_dim <- as.integer(params$m_manifold_dim)
  }
  
  if (params$max_iter <= 0) {
    stop("max_iter must be positive")
  }
  
  if (params$tolerance <= 0) {
    stop("tolerance must be positive")
  }
  
  if (params$lambda_manifold < 0) {
    stop("lambda_manifold must be non-negative")
  }
  
  if (params$lambda_spatial_smooth < 0) {
    stop("lambda_spatial_smooth must be non-negative")
  }
  
  if (params$lambda_gamma < 0) {
    stop("lambda_gamma must be non-negative")
  }
  
  if (params$num_neighbors_mfd < 1) {
    stop("num_neighbors_mfd must be at least 1")
  }
  
  if (params$num_neighbors_Lsp < 1) {
    stop("num_neighbors_Lsp must be at least 1")
  }
  
  if (params$n_cores < 1) {
    stop("n_cores must be at least 1")
  }
  
  if (!(params$verbose_level %in% c(0, 1, 2))) {
    stop("verbose_level must be 0, 1, or 2")
  }
  
  structure(params, class = "manifold_control", preset = preset)
}

#' Print method for manifold_control objects
#'
#' @param x A manifold_control object
#' @param ... Ignored
#' @export
print.manifold_control <- function(x, ...) {
  preset <- attr(x, "preset")
  cat("Manifold HRF Control Parameters (preset: ", preset, ")\n", sep = "")
  cat("----------------------------------------\n")
  
  # Group parameters by category
  cat("Manifold construction:\n")
  cat("  m_manifold_dim:", ifelse(is.null(x$m_manifold_dim), "auto", x$m_manifold_dim), "\n")
  cat("  lambda_manifold:", x$lambda_manifold, "\n")
  cat("  num_neighbors_mfd:", x$num_neighbors_mfd, "\n")
  
  cat("\nSpatial smoothing:\n")
  cat("  lambda_spatial_smooth:", x$lambda_spatial_smooth, "\n")
  cat("  num_neighbors_Lsp:", x$num_neighbors_Lsp, "\n")
  
  cat("\nOptimization:\n")
  cat("  max_iter:", x$max_iter, "\n")
  cat("  tolerance:", x$tolerance, "\n")
  cat("  lambda_gamma:", x$lambda_gamma, "\n")
  cat("  lambda_condition_betas:", x$lambda_condition_betas, "\n")
  cat("  orthogonal_approx:", x$orthogonal_approx, "\n")
  
  cat("\nComputational:\n")
  cat("  n_cores:", x$n_cores, "\n")
  cat("  verbose_level:", x$verbose_level, "\n")
  cat("  ram_threshold_GB:", x$ram_threshold_GB, "\n")
  cat("  generate_qc_report:", x$generate_qc_report, "\n")
  
  invisible(x)
}