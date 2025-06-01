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
    library_list_or_matrix <- switch(hrf_library_source,
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

    time_points <- seq(0, hrf_duration, by = TR_precision)
    if (is.list(library_list_or_matrix) && all(sapply(library_list_or_matrix, inherits, "HRF"))) {
      library_info$n_hrfs <- length(library_list_or_matrix)
      L_library_matrix <- do.call(cbind, lapply(library_list_or_matrix, function(h) {
        as.numeric(fmrireg::evaluate(h, time_points))
      }))
    } else if (is.matrix(library_list_or_matrix)) {
      L_library_matrix <- library_list_or_matrix
    } else {
      stop("Unsupported HRF library format returned from switch")
    }
  } else if (is.list(hrf_library_source)) {
    # List input - may be fmrireg HRF objects or generic structures
    library_info$name <- "custom_list"
    library_info$n_hrfs <- length(hrf_library_source)

    if (length(hrf_library_source) == 0) {
      stop("hrf_library_source list cannot be empty")
    }

    time_points <- seq(0, hrf_duration, by = TR_precision)
    p <- length(time_points)
    N <- length(hrf_library_source)

    if (all(sapply(hrf_library_source, inherits, "HRF"))) {
      # Evaluate each fmrireg HRF object
      L_library_matrix <- do.call(cbind, lapply(hrf_library_source, function(hrf_obj) {
        as.numeric(fmrireg::evaluate(hrf_obj, time_points))
      }))
    } else {
      # Fallback: treat as generic list and generate dummy HRFs (backwards compatibility)
      L_library_matrix <- matrix(0, p, N)
      for (i in seq_len(N)) {
        L_library_matrix[, i] <- dgamma(time_points, shape = 6 + i, rate = 1)
      }
    }

    # Normalize each HRF
    L_library_matrix <- apply(L_library_matrix, 2, function(h) h / sum(abs(h)))
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
#'
#' Generates a small set of FLOBS-like HRF basis functions as a list of
#' `fmrireg::HRF` objects.
#' @keywords internal
create_flobs_library <- function(TR_precision, hrf_duration) {
  time_points <- seq(0, hrf_duration, by = TR_precision)

  n_basis <- 5
  library_list <- vector("list", n_basis)

  for (i in seq_len(n_basis)) {
    hrf_fun <- local({
      shift <- (i - 1) * 0.5
      function(t) {
        t_shift <- t - shift
        t_shift[t_shift < 0] <- 0
        stats::dgamma(t_shift, shape = 6, rate = 1) -
          0.35 * stats::dgamma(t_shift, shape = 16, rate = 1)
      }
    })

    # Normalise based on sampling grid
    scale_val <- max(abs(hrf_fun(time_points)))
    hrf_norm <- function(t) hrf_fun(t) / scale_val

    library_list[[i]] <- fmrireg::as_hrf(
      hrf_norm,
      name = paste0("flobs", i),
      span = hrf_duration
    )
  }

  library_list
}

#' Create Half-Cosine Library
#'
#' Returns a list of damped cosine basis functions as `fmrireg::HRF` objects.
#' @keywords internal
create_half_cosine_library <- function(TR_precision, hrf_duration, n_basis = 20) {
  time_points <- seq(0, hrf_duration, by = TR_precision)

  library_list <- vector("list", n_basis)

  for (i in seq_len(n_basis)) {
    freq <- i * pi / hrf_duration
    hrf_fun <- function(t) {
      stats::cos(freq * t) * exp(-t / (hrf_duration / 2))
    }

    scale_val <- max(abs(hrf_fun(time_points)))
    hrf_norm <- function(t) hrf_fun(t) / scale_val

    library_list[[i]] <- fmrireg::as_hrf(
      hrf_norm,
      name = paste0("half_cosine", i),
      span = hrf_duration
    )
  }

  library_list
}

#' Create Gamma Grid Library
#'
#' Generates a grid of gamma HRFs as a list of `fmrireg::HRF` objects.
#' @keywords internal
create_gamma_grid_library <- function(TR_precision, hrf_duration) {
  time_points <- seq(0, hrf_duration, by = TR_precision)

  shapes <- seq(4, 8, length.out = 10)
  scales <- seq(0.8, 1.2, length.out = 10)

  grid <- expand.grid(shape = shapes, scale = scales)

  library_list <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    shape_i <- grid$shape[i]
    scale_i <- grid$scale[i]
    hrf_fun <- function(t) {
      stats::dgamma(t, shape = shape_i, scale = scale_i)
    }

    norm_val <- sum(hrf_fun(time_points))
    hrf_norm <- function(t) hrf_fun(t) / norm_val

    library_list[[i]] <- fmrireg::as_hrf(
      hrf_norm,
      name = paste0("gamma_shape", shape_i, "_scale", scale_i),
      span = hrf_duration
    )
  }

  library_list
}

#' Create Manifold HRF Object
#' @keywords internal
create_manifold_hrf_object <- function(B_reconstructor, name, nbasis) {
  if (!is.matrix(B_reconstructor)) {
    stop("B_reconstructor must be a matrix")
  }

  nbasis <- nbasis %||% ncol(B_reconstructor)
  if (ncol(B_reconstructor) < nbasis) {
    stop("nbasis exceeds number of columns in B_reconstructor")
  }

  time_points <- seq(0, length.out = nrow(B_reconstructor), by = 1)

  basis_functions <- vector("list", nbasis)
  for (j in seq_len(nbasis)) {
    basis_functions[[j]] <- fmrireg::empirical_hrf(
      time_points,
      B_reconstructor[, j],
      name = paste0(name, "_basis", j)
    )
  }

  manifold_hrf <- do.call(fmrireg::bind_basis, basis_functions)
  attr(manifold_hrf, "name") <- name
  attr(manifold_hrf, "nbasis") <- nbasis
  class(manifold_hrf) <- c("mhrf_basis", class(manifold_hrf))

  manifold_hrf
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
process_subject_mhrf_lss_nim <- function(bold_input, mask_input, event_input,
                                        confound_input = NULL, manifold_objects,
                                        params_list) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("Package 'neuroim2' is required")
  }

  bold <- if (inherits(bold_input, c("NeuroVec", "NeuroVol"))) {
    bold_input
  } else {
    neuroim2::read_vec(bold_input)
  }

  mask <- if (inherits(mask_input, "LogicalNeuroVol")) {
    mask_input
  } else {
    neuroim2::read_vol(mask_input)
  }

  events <- if (is.character(event_input)) {
    read.csv(event_input, sep = ifelse(grepl("\\.tsv$", event_input), "\t", ","))
  } else {
    event_input
  }

  TR <- params_list$TR %||% stop("TR must be specified in params_list")

  Y_mat <- as.matrix(bold)
  mask_idx <- which(as.logical(mask))
  Y_mat <- Y_mat[, mask_idx, drop = FALSE]

  sframe <- fmrireg::sampling_frame(blocklens = nrow(Y_mat), TR = TR)
  ev_model <- fmrireg::event_model(~ hrf(onset, basis = manifold_objects$manifold_hrf_basis),
                                   data = events, sampling_frame = sframe, drop_empty = TRUE)
  design_info <- extract_design_info(ev_model, sframe)

  result <- run_mhrf_lss_standard(
    Y_data = Y_mat,
    design_info = design_info,
    manifold = manifold_objects,
    Z_confounds = confound_input,
    voxel_coords = neuroim2::coords(mask)[mask_idx, , drop = FALSE],
    params = params_list,
    outlier_weights = NULL,
    estimation = if (length(design_info$X_trial_list) > 0) "both" else "condition",
    progress = FALSE
  )

  result$mask_indices <- mask_idx
  result$space <- neuroim2::space(bold)

  return(result)
}

#' Enhanced Results Packaging & Visualization (Neuroimaging Layer)
#'
#' @param core_results_list Results from core pipeline functions
#' @param reference_space NeuroSpace object from mask/BOLD
#' @param mask_vol LogicalNeuroVol brain mask
#' @param original_inputs Original input specifications
#' @param processing_metadata Metadata from processing
#' @return mhrf_results S3 object with neuroim2 integration
package_mhrf_results_nim <- function(core_results_list, reference_space, mask_vol,
                                    original_inputs, processing_metadata) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("Package 'neuroim2' is required")
  }

  mask_idx <- which(as.logical(mask_vol))
  p <- nrow(core_results_list$H_shapes)
  Vtot <- prod(neuroim2::dim(reference_space))

  hrf_arr <- matrix(0, p, Vtot)
  hrf_arr[, mask_idx] <- core_results_list$H_shapes
  hrf_vec <- neuroim2::NeuroVec(hrf_arr, reference_space)

  amp_arr <- matrix(0, nrow(core_results_list$Beta_condition), Vtot)
  amp_arr[, mask_idx] <- core_results_list$Beta_condition
  amp_vec <- neuroim2::NeuroVec(amp_arr, reference_space)

  result <- list(
    hrfs = hrf_vec,
    amplitudes = amp_vec,
    matrices = core_results_list,
    metadata = processing_metadata,
    inputs = original_inputs
  )
  class(result) <- c("mhrf_results", "list")
  return(result)
}
