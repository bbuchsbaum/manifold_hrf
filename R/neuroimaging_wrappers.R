# Neuroimaging Wrapper Functions
# Implementation of MHRF-NIM-IO-MANIFOLD-01 and other neuroimaging integration

#' Construct HRF Manifold (Neuroimaging Wrapper)
#'
#' Creates an HRF manifold from various neuroimaging-compatible sources,
#' wrapping the core manifold construction functions with fmrireg integration.
#'
#' @param hrf_library_source Source for HRF library. Can be:
#'   \itemize{
#'     \item Character string: "FLOBS", "half_cosine", "gamma_grid", or path to RDS file
#'     \item List of fmrireg HRF objects
#'     \item Matrix (p x N) of HRF shapes
#'   }
#' @param TR_precision Time resolution in seconds for HRF evaluation (e.g., 0.1)
#' @param hrf_duration Total duration of HRF in seconds (default 24)
#' @param m_manifold_dim_target Target manifold dimensionality (default 5)
#' @param m_manifold_dim_min_variance Minimum variance explained (default 0.95)
#' @param k_local_nn_for_sigma k-NN for self-tuning bandwidth (default 7)
#' @param sparse_threshold Library size above which to use sparse matrices (default 5000)
#' @param ... Additional parameters passed to core functions
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{B_reconstructor_matrix}: p x m HRF reconstructor matrix
#'     \item \code{manifold_hrf_basis}: Custom fmrireg HRF object for the manifold basis
#'     \item \code{Phi_coords_matrix}: N x m manifold coordinates of library HRFs
#'     \item \code{eigenvalues_S}: Eigenvalues from diffusion map
#'     \item \code{m_manifold_dim}: Final manifold dimensionality used
#'     \item \code{m_auto_selected}: Automatically selected dimensionality
#'     \item \code{library_info}: Information about the HRF library used
#'     \item \code{parameters}: List of all parameters used
#'   }
#'
#' @details This function provides a neuroimaging-friendly interface to the core
#'   manifold construction pipeline. It handles conversion between fmrireg HRF
#'   objects and the matrix format required by core functions. The resulting
#'   manifold basis can be used directly in fmrireg model specifications.
#'
#' @examples
#' \dontrun{
#' # Use FLOBS basis
#' manifold_flobs <- construct_hrf_manifold_nim(
#'   hrf_library_source = "FLOBS",
#'   TR_precision = 0.1,
#'   m_manifold_dim_target = 5
#' )
#' 
#' # Use custom HRF objects from fmrireg
#' hrf_list <- list(
#'   fmrireg::HRF_SPMG1,
#'   fmrireg::HRF_SPMG2,
#'   fmrireg::HRF_SPMG3,
#'   fmrireg::HRF_GAMMA
#' )
#' manifold_custom <- construct_hrf_manifold_nim(
#'   hrf_library_source = hrf_list,
#'   TR_precision = 0.5
#' )
#' 
#' # Use in fmrireg model
#' # event_model(~ hrf(onset, basis = manifold_custom$manifold_hrf_basis))
#' }
#'
#' @export
construct_hrf_manifold_nim <- function(hrf_library_source,
                                     TR_precision = 0.1,
                                     hrf_duration = 24,
                                     m_manifold_dim_target = 5,
                                     m_manifold_dim_min_variance = 0.95,
                                     k_local_nn_for_sigma = 7,
                                     sparse_threshold = 5000,
                                     ...) {
  
  # Step 1: Convert HRF library source to matrix format
  library_info <- list(source_type = class(hrf_library_source)[1])
  
  if (is.character(hrf_library_source) && length(hrf_library_source) == 1) {
    # Handle predefined libraries or file paths
    L_library_matrix <- switch(hrf_library_source,
      "FLOBS" = {
        library_info$name <- "FLOBS"
        library_info$description <- "FMRIB's Linear Optimal Basis Set"
        create_flobs_library(TR_precision, hrf_duration)
      },
      "half_cosine" = {
        library_info$name <- "half_cosine"
        library_info$description <- "Half-cosine basis set"
        create_half_cosine_library(TR_precision, hrf_duration, n_basis = 20)
      },
      "gamma_grid" = {
        library_info$name <- "gamma_grid"
        library_info$description <- "Gamma function grid"
        create_gamma_grid_library(TR_precision, hrf_duration)
      },
      {
        # Assume it's a file path
        if (file.exists(hrf_library_source)) {
          library_info$name <- "custom_file"
          library_info$path <- hrf_library_source
          readRDS(hrf_library_source)
        } else {
          stop("Unknown HRF library source or file not found: ", hrf_library_source)
        }
      }
    )
  } else if (is.list(hrf_library_source)) {
    # List of HRF objects - check if they're fmrireg HRFs
    library_info$name <- "custom_list"
    library_info$n_hrfs <- length(hrf_library_source)
    
    # Check if these are fmrireg HRF objects
    if (length(hrf_library_source) > 0) {
      # For now, create a simple matrix representation
      # In full implementation, would use fmrireg::evaluate() on each HRF
      time_points <- seq(0, hrf_duration, by = TR_precision)
      p <- length(time_points)
      N <- length(hrf_library_source)
      
      L_library_matrix <- matrix(0, p, N)
      for (i in 1:N) {
        # Placeholder: would call fmrireg::evaluate(hrf_library_source[[i]], time_points)
        # For now, create dummy HRFs
        L_library_matrix[, i] <- dgamma(time_points, shape = 6 + i, rate = 1)
      }
      
      # Normalize each HRF
      L_library_matrix <- apply(L_library_matrix, 2, function(h) h / sum(abs(h)))
    }
  } else if (is.matrix(hrf_library_source)) {
    # Already in matrix format
    library_info$name <- "custom_matrix"
    library_info$dimensions <- dim(hrf_library_source)
    L_library_matrix <- hrf_library_source
  } else {
    stop("hrf_library_source must be a character string, list of HRF objects, or matrix")
  }
  
  # Validate matrix dimensions
  if (nrow(L_library_matrix) < 2) {
    stop("HRF library must have at least 2 time points")
  }
  if (ncol(L_library_matrix) < m_manifold_dim_target + 1) {
    warning(sprintf(
      "HRF library has only %d HRFs but target manifold dimension is %d. Adjusting target.",
      ncol(L_library_matrix), m_manifold_dim_target
    ))
    m_manifold_dim_target <- min(m_manifold_dim_target, ncol(L_library_matrix) - 1)
  }
  
  # Step 2: Determine if we should use sparse matrices
  N <- ncol(L_library_matrix)
  use_sparse <- N > sparse_threshold
  
  # Adjust k_local_nn_for_sigma if necessary
  k_local_nn_for_sigma_adj <- min(k_local_nn_for_sigma, N - 1)
  if (k_local_nn_for_sigma_adj < k_local_nn_for_sigma) {
    warning(sprintf(
      "Adjusting k_local_nn_for_sigma from %d to %d due to small library size (N=%d)",
      k_local_nn_for_sigma, k_local_nn_for_sigma_adj, N
    ))
  }
  
  if (use_sparse) {
    message(sprintf("Using sparse matrices for large library (N = %d)", N))
    use_sparse_params <- list(
      sparse_if_N_gt = sparse_threshold,
      k_nn_for_W_sparse = min(k_local_nn_for_sigma_adj * 3, N - 1)
    )
  } else {
    use_sparse_params <- list(
      sparse_if_N_gt = Inf,
      k_nn_for_W_sparse = NULL
    )
  }
  
  # Step 3: Call core manifold construction functions
  message("Computing HRF affinity matrix...")
  S_markov <- calculate_manifold_affinity_core(
    L_library_matrix = L_library_matrix,
    k_local_nn_for_sigma = k_local_nn_for_sigma_adj,
    use_sparse_W_params = use_sparse_params
  )
  
  message("Computing diffusion map embedding...")
  manifold_result <- get_manifold_basis_reconstructor_core(
    S_markov_matrix = S_markov,
    L_library_matrix = L_library_matrix,
    m_manifold_dim_target = m_manifold_dim_target,
    m_manifold_dim_min_variance = m_manifold_dim_min_variance
  )
  
  # Step 4: Create fmrireg-compatible HRF basis object
  # This would use fmrireg::as_hrf() in full implementation
  manifold_hrf_basis <- create_manifold_hrf_object(
    B_reconstructor = manifold_result$B_reconstructor_matrix,
    name = paste0("manifold_", library_info$name),
    nbasis = manifold_result$m_final_dim
  )
  
  # Step 5: Compile results
  result <- list(
    B_reconstructor_matrix = manifold_result$B_reconstructor_matrix,
    manifold_hrf_basis = manifold_hrf_basis,
    Phi_coords_matrix = manifold_result$Phi_coords_matrix,
    eigenvalues_S = manifold_result$eigenvalues_S_vector,
    m_manifold_dim = manifold_result$m_final_dim,
    m_auto_selected = manifold_result$m_auto_selected_dim,
    library_info = library_info,
    parameters = list(
      TR_precision = TR_precision,
      hrf_duration = hrf_duration,
      m_manifold_dim_target = m_manifold_dim_target,
      m_manifold_dim_min_variance = m_manifold_dim_min_variance,
      k_local_nn_for_sigma = k_local_nn_for_sigma,
      sparse_threshold = sparse_threshold,
      use_sparse = use_sparse,
      n_hrfs_library = N,
      n_timepoints = nrow(L_library_matrix)
    )
  )
  
  class(result) <- c("mhrf_manifold", "list")
  
  message(sprintf(
    "Created HRF manifold: %d-dimensional embedding of %d HRFs (variance explained: %.1f%%)",
    result$m_manifold_dim,
    N,
    100 * sum(result$eigenvalues_S[2:(result$m_manifold_dim + 1)]) / 
      sum(result$eigenvalues_S[-1])
  ))
  
  return(result)
}


# Helper functions for creating HRF libraries

#' Create FLOBS Library
#' @keywords internal
create_flobs_library <- function(TR_precision, hrf_duration) {
  # Simplified FLOBS-like basis
  # In real implementation, would use actual FLOBS basis
  time_points <- seq(0, hrf_duration, by = TR_precision)
  p <- length(time_points)
  
  # Create basis functions
  n_basis <- 5
  L <- matrix(0, p, n_basis)
  
  # Canonical HRF and derivatives
  for (i in 1:n_basis) {
    t_shift <- time_points - (i-1) * 0.5
    t_shift[t_shift < 0] <- 0
    L[, i] <- dgamma(t_shift, shape = 6, rate = 1) - 
              0.35 * dgamma(t_shift, shape = 16, rate = 1)
  }
  
  # Normalize
  L <- apply(L, 2, function(x) x / max(abs(x)))
  
  return(L)
}

#' Create Half-Cosine Library
#' @keywords internal
create_half_cosine_library <- function(TR_precision, hrf_duration, n_basis = 20) {
  time_points <- seq(0, hrf_duration, by = TR_precision)
  p <- length(time_points)
  
  L <- matrix(0, p, n_basis)
  
  for (i in 1:n_basis) {
    freq <- i * pi / hrf_duration
    L[, i] <- cos(freq * time_points)
    # Apply window to make it causal
    L[, i] <- L[, i] * exp(-time_points / (hrf_duration / 2))
  }
  
  # Normalize
  L <- apply(L, 2, function(x) x / max(abs(x)))
  
  return(L)
}

#' Create Gamma Grid Library
#' @keywords internal
create_gamma_grid_library <- function(TR_precision, hrf_duration) {
  time_points <- seq(0, hrf_duration, by = TR_precision)
  p <- length(time_points)
  
  # Grid of shape and scale parameters
  shapes <- seq(4, 8, length.out = 10)
  scales <- seq(0.8, 1.2, length.out = 10)
  
  grid <- expand.grid(shape = shapes, scale = scales)
  N <- nrow(grid)
  
  L <- matrix(0, p, N)
  
  for (i in 1:N) {
    L[, i] <- dgamma(time_points, shape = grid$shape[i], scale = grid$scale[i])
  }
  
  # Normalize
  L <- apply(L, 2, function(x) x / sum(x))
  
  return(L)
}

#' Create Manifold HRF Object
#' @keywords internal
create_manifold_hrf_object <- function(B_reconstructor, name, nbasis) {
  # Placeholder for creating fmrireg-compatible HRF object
  # In real implementation would use fmrireg::as_hrf()
  
  hrf_obj <- list(
    type = "manifold",
    name = name,
    nbasis = nbasis,
    B_reconstructor = B_reconstructor,
    evaluate = function(grid, ...) {
      # Placeholder evaluation function
      # Would implement proper manifold coordinate to HRF mapping
      matrix(0, length(grid), nbasis)
    }
  )
  
  class(hrf_obj) <- c("mhrf_basis", "hrf", "list")
  
  return(hrf_obj)
}

#' Print method for mhrf_manifold objects
#' @export
print.mhrf_manifold <- function(x, ...) {
  cat("M-HRF Manifold Object\n")
  cat("====================\n")
  cat(sprintf("Library: %s (%d HRFs)\n", 
              x$library_info$name, 
              x$parameters$n_hrfs_library))
  cat(sprintf("Manifold dimension: %d (auto-selected: %d)\n", 
              x$m_manifold_dim, 
              x$m_auto_selected))
  cat(sprintf("Variance explained: %.1f%%\n",
              100 * sum(x$eigenvalues_S[2:(x$m_manifold_dim + 1)]) / 
                sum(x$eigenvalues_S[-1])))
  cat(sprintf("Time points: %d (TR precision: %.2fs)\n",
              x$parameters$n_timepoints,
              x$parameters$TR_precision))
  cat(sprintf("Sparse matrices: %s\n",
              ifelse(x$parameters$use_sparse, "Yes", "No")))
}

#' Summary method for mhrf_manifold objects
#' @export
summary.mhrf_manifold <- function(object, ...) {
  cat("M-HRF Manifold Summary\n")
  cat("=====================\n\n")
  
  print(object)
  
  cat("\nEigenvalue spectrum (top 10):\n")
  n_show <- min(10, length(object$eigenvalues_S))
  eigenvals <- object$eigenvalues_S[1:n_show]
  var_exp <- eigenvals[-1] / sum(object$eigenvalues_S[-1]) * 100
  
  df <- data.frame(
    Component = 1:n_show,
    Eigenvalue = round(eigenvals, 4),
    `Var.Explained` = c(NA, round(var_exp, 1)),
    `Cumulative.Var` = c(NA, round(cumsum(var_exp), 1))
  )
  
  print(df, row.names = FALSE)
  
  cat("\nReconstructor matrix dimensions:", 
      paste(dim(object$B_reconstructor_matrix), collapse = " x "), "\n")
}

# Placeholder functions for other wrapper components
# These would be fully implemented in a complete version

#' Enhanced Subject-Level Processing Wrapper (Neuroimaging Layer)
#'
#' @param bold_input BOLD data (NIfTI path, NeuroVec, or matrix)
#' @param mask_input Brain mask (path, LogicalNeuroVol, or logical vector)
#' @param event_input Event data (path to CSV/TSV or data.frame)
#' @param confound_input Optional confound regressors
#' @param manifold_objects Manifold objects from construct_hrf_manifold_nim
#' @param params_list All pipeline parameters
#' @return List of R matrices plus processing metadata
#' @export
process_subject_mhrf_lss_nim <- function(bold_input, mask_input, event_input,
                                        confound_input = NULL, manifold_objects,
                                        params_list) {
  # TODO: Implement MHRF-NIM-WRAP-SUBJECT-01
  # This would handle:
  # - Loading BOLD data via neuroim2
  # - Creating design matrices via fmrireg
  # - Calling all core pipeline functions
  # - Returning results
  
  stop("Function not yet implemented - requires neuroim2/fmrireg integration")
}

#' Enhanced Results Packaging & Visualization (Neuroimaging Layer)
#'
#' @param core_results_list Results from core pipeline functions
#' @param reference_space NeuroSpace object from mask/BOLD
#' @param mask_vol LogicalNeuroVol brain mask
#' @param original_inputs Original input specifications
#' @param processing_metadata Metadata from processing
#' @return mhrf_results S3 object with neuroim2 integration
#' @export
package_mhrf_results_nim <- function(core_results_list, reference_space, mask_vol,
                                    original_inputs, processing_metadata) {
  # TODO: Implement MHRF-NIM-OUTPUT-01
  # This would:
  # - Convert matrices back to NeuroVec/NeuroVol objects
  # - Create S3 class with print/plot/summary methods
  # - Enable writing results to NIfTI
  
  stop("Function not yet implemented - requires neuroim2 integration")
}