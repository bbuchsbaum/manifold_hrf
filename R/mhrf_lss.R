#' Analyze fMRI Data with M-HRF-LSS
#'
#' Main unified function for running the complete M-HRF-LSS pipeline. This function
#' provides a simple interface to estimate voxel-specific HRFs and trial-wise
#' amplitudes from fMRI data.
#'
#' @param Y_data An n x V numeric matrix of fMRI time series (n timepoints, V voxels),
#'   or a NeuroVec object from the neuroim2 package
#' @param events A data frame with columns:
#'   \itemize{
#'     \item \code{onset}: Event onset times in seconds
#'     \item \code{duration}: Event duration in seconds (default: 0)
#'     \item \code{condition}: Factor or character vector of condition labels
#'     \item \code{trial_type}: Optional trial identifiers for trial-wise estimation
#'   }
#' @param TR Repetition time in seconds (default: 2)
#' @param preset Character string specifying parameter preset: "conservative",
#'   "balanced" (default), "aggressive", "fast", "quality", or "robust"
#' @param hrf_library Source for HRF library. Options:
#'   \itemize{
#'     \item \code{"auto"}: Automatically generate appropriate library (default)
#'     \item \code{"spmg1"}: SPM canonical HRF variants
#'     \item \code{"gamma"}: Gamma function variants  
#'     \item \code{"custom"}: User-provided p x N matrix
#'   }
#' @param voxel_mask Optional logical vector or 3D array indicating which voxels
#'   to analyze. If NULL, all non-zero variance voxels are analyzed.
#' @param n_jobs Number of parallel jobs for voxel processing (default: 1).
#'   Set to -1 to use all available cores.
#' @param verbose Logical (default: TRUE) or integer (0-3) controlling output:
#'   \itemize{
#'     \item \code{0}: Silent
#'     \item \code{1}/\code{TRUE}: Progress milestones
#'     \item \code{2}: Detailed progress
#'     \item \code{3}: Debug output
#'   }
#' @param save_intermediate Logical; whether to save intermediate results for
#'   debugging (default: FALSE)
#' @param output_dir Directory for saving results and intermediate files
#'   (default: temporary directory)
#' @param ... Additional parameters to override preset values. See 
#'   \code{\link{get_preset_params}} for available options.
#'
#' @return An object of class "mhrf_result" containing:
#'   \itemize{
#'     \item \code{hrf_shapes}: p x V matrix of estimated HRF shapes
#'     \item \code{amplitudes}: k x V matrix of condition-level amplitudes
#'     \item \code{trial_amplitudes}: List of trial-wise amplitude estimates
#'     \item \code{manifold_coords}: m x V matrix of manifold coordinates
#'     \item \code{qc_metrics}: Quality control metrics and diagnostics
#'     \item \code{metadata}: Processing details and parameters used
#'     \item \code{call}: The original function call
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic usage with data frame events
#' result <- mhrf_lss(
#'   Y_data = fmri_matrix,
#'   events = event_df,
#'   TR = 2
#' )
#' 
#' # Use robust preset for noisy data
#' result <- mhrf_lss(
#'   Y_data = fmri_matrix,
#'   events = event_df,
#'   TR = 2,
#'   preset = "robust"
#' )
#' 
#' # Parallel processing with custom parameters
#' result <- mhrf_lss(
#'   Y_data = fmri_matrix,
#'   events = event_df,
#'   TR = 2,
#'   n_jobs = 4,
#'   lambda_gamma = 0.05,
#'   m_manifold_dim_target = 6
#' )
#' 
#' # With brain mask
#' result <- mhrf_lss(
#'   Y_data = fmri_matrix,
#'   events = event_df,
#'   TR = 2,
#'   voxel_mask = brain_mask
#' )
#' }
#'
#' @export
#' @seealso 
#' \code{\link{get_preset_params}} for parameter presets,
#' \code{\link{summary.mhrf_result}} for result summary,
#' \code{\link{plot.mhrf_result}} for diagnostic plots
mhrf_analyze <- function(Y_data,
                     events,
                     TR = 2,
                     preset = "balanced",
                     hrf_library = "auto",
                     voxel_mask = NULL,
                     n_jobs = 1,
                     verbose = TRUE,
                     save_intermediate = FALSE,
                     output_dir = tempdir(),
                     ...) {
  
  # Capture call for reproducibility
  mc <- match.call()
  start_time <- Sys.time()
  
  # Set up verbosity
  verbose_level <- if (is.logical(verbose)) {
    if (verbose) 1L else 0L
  } else {
    as.integer(verbose)
  }
  
  # Create progress tracker
  progress <- .create_progress_tracker(verbose_level)
  progress$start("M-HRF-LSS Analysis")
  
  # Step 1: Input validation and preprocessing
  progress$update("Validating inputs...")
  
  # Comprehensive input validation
  data_validation <- .validate_Y_data(Y_data)
  events_validation <- .validate_events(events, data_validation$n_timepoints, TR)
  .validate_parameters(TR, preset, data_validation$n_voxels, list(...))
  mask_validation <- .validate_voxel_mask(voxel_mask, data_validation$n_voxels)
  .check_system_requirements(data_validation$n_voxels, data_validation$n_timepoints, preset)
  
  # Convert inputs to standard format
  data_info <- .prepare_data_inputs(
    Y_data = Y_data,
    events = events_validation$events_validated,
    TR = TR,
    voxel_mask = voxel_mask,
    progress = progress
  )
  
  Y_matrix <- data_info$Y_matrix
  n_timepoints <- data_info$n_timepoints
  n_voxels <- data_info$n_voxels
  voxel_indices <- data_info$voxel_indices
  
  # Step 2: Load parameters
  progress$update("Loading parameters...")
  
  # Get base parameters from preset
  params <- get_preset_params(preset, n_voxels = n_voxels)
  
  # Override with user parameters
  user_params <- list(...)
  for (param_name in names(user_params)) {
    params[[param_name]] <- user_params[[param_name]]
  }
  
  # Add system parameters
  params$TR <- TR
  params$n_jobs <- .determine_n_jobs(n_jobs)
  params$verbose_level <- verbose_level
  params$save_intermediate <- save_intermediate
  params$output_dir <- output_dir
  
  # Add default values for parameters that might be missing
  params$p_hrf <- params$p_hrf %||% 25
  params$screen_voxels <- params$screen_voxels %||% TRUE
  params$apply_hrf_constraints <- params$apply_hrf_constraints %||% FALSE
  params$hrf_peak_range <- params$hrf_peak_range %||% c(2, 10)
  params$estimate_trials <- params$estimate_trials %||% TRUE
  params$generate_qc_report <- params$generate_qc_report %||% FALSE
  params$use_robust_manifold <- params$use_robust_manifold %||% params$use_robust_svd
  params$fallback_to_pca <- params$fallback_to_pca %||% FALSE
  params$adaptive_smoothing <- params$adaptive_smoothing %||% FALSE
  params$edge_preserve <- params$edge_preserve %||% FALSE
  params$chunk_size <- params$chunk_size %||% 1000
  params$orthogonal_approx <- params$orthogonal_approx %||% FALSE
  params$data_checks <- params$data_checks %||% TRUE
  
  # Step 3: Create design matrices
  progress$update("Creating design matrices...")
  
  design_info <- .create_design_matrices(
    events = events,
    n_timepoints = n_timepoints,
    TR = TR,
    hrf_length = params$p_hrf
  )
  
  X_condition_list <- design_info$X_condition_list
  X_trial_list <- design_info$X_trial_list
  n_conditions <- design_info$n_conditions
  n_trials <- design_info$n_trials
  
  # Validate parameters (now that we have design matrices)
  if (params$data_checks) {
    params$validate_data(Y_matrix, X_condition_list)
  }
  
  # Step 4: Create or load HRF library
  progress$update("Creating HRF library...")
  
  hrf_lib_info <- .prepare_hrf_library(
    hrf_library = hrf_library,
    p_hrf = params$p_hrf,
    TR = TR,
    params = params
  )
  
  L_library <- hrf_lib_info$L_library
  n_hrfs <- hrf_lib_info$n_hrfs
  
  # Step 5: Handle problematic voxels
  if (params$screen_voxels) {
    progress$update("Screening voxels...")
    
    screening_result <- screen_voxels(Y_matrix)
    valid_voxels <- which(screening_result$keep)
    
    if (length(valid_voxels) < n_voxels) {
      progress$message(sprintf("Excluding %d problematic voxels", 
                              n_voxels - length(valid_voxels)))
      Y_matrix <- Y_matrix[, valid_voxels, drop = FALSE]
      voxel_indices <- voxel_indices[valid_voxels]
      n_voxels <- length(valid_voxels)
    }
  }
  
  # Step 6: Run core M-HRF-LSS pipeline
  progress$update("Running M-HRF-LSS pipeline...")
  
  # Component 0: HRF Manifold Construction
  progress$update("  Component 0: Constructing HRF manifold...", level = 2)
  
  manifold_result <- .run_manifold_construction(
    L_library = L_library,
    params = params,
    progress = progress
  )
  
  # Component 1: Voxel-wise HRF Estimation
  progress$update("  Component 1: Estimating voxel-wise HRFs...", level = 2)
  
  voxelwise_result <- .run_voxelwise_estimation(
    Y_data = Y_matrix,
    X_condition_list = X_condition_list,
    manifold = manifold_result,
    params = params,
    progress = progress
  )
  
  # Component 2: Spatial Smoothing
  progress$update("  Component 2: Applying spatial smoothing...", level = 2)
  
  smoothing_result <- .run_spatial_smoothing(
    Xi_matrix = voxelwise_result$Xi_ident,
    voxel_coords = data_info$voxel_coords,
    params = params,
    progress = progress,
    Y_data = Y_matrix
  )
  
  # Component 3: HRF Reconstruction
  progress$update("  Component 3: Reconstructing HRF shapes...", level = 2)
  
  hrf_shapes <- reconstruct_hrf_shapes_core(
    B_reconstructor_matrix = manifold_result$B_reconstructor,
    Xi_smoothed_matrix = smoothing_result$Xi_smooth
  )
  
  # Apply physiological constraints if requested
  if (params$apply_hrf_constraints) {
    hrf_constraint_result <- apply_hrf_physiological_constraints(
      hrf_matrix = hrf_shapes,
      TR = TR,
      peak_range = params$hrf_peak_range,
      enforce_positive = TRUE
    )
    hrf_shapes <- hrf_constraint_result$hrf_constrained
  }
  
  # Component 4: Trial-wise LSS (if requested)
  trial_amplitudes <- NULL
  if (n_trials > 0 && params$estimate_trials) {
    progress$update("  Component 4: Estimating trial-wise amplitudes...", level = 2)
    
    trial_amplitudes <- .run_trial_estimation(
      Y_data = Y_matrix,
      X_trial_list = X_trial_list,
      hrf_shapes = hrf_shapes,
      params = params,
      progress = progress
    )
  }
  
  # Step 7: Compute QC metrics
  progress$update("Computing quality metrics...")
  
  qc_metrics <- .compute_qc_metrics(
    Y_data = Y_matrix,
    X_condition_list = X_condition_list,
    hrf_shapes = hrf_shapes,
    amplitudes = voxelwise_result$Beta_ident,
    manifold_coords = smoothing_result$Xi_smooth,
    params = params
  )
  
  # Step 8: Package results
  progress$update("Packaging results...")
  
  # Create metadata
  metadata <- list(
    n_timepoints = n_timepoints,
    n_voxels_input = data_info$n_voxels_original,
    n_voxels_analyzed = n_voxels,
    n_conditions = n_conditions,
    n_trials = n_trials,
    n_hrfs_library = n_hrfs,
    manifold_dim = manifold_result$m_final,
    manifold_method = manifold_result$method_used,
    preset_used = preset,
    parameters = params,
    runtime_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
    version = packageVersion("manifoldhrf")
  )
  
  # Add condition names to amplitudes matrix
  amplitudes_matrix <- voxelwise_result$Beta_ident
  if (is.null(rownames(amplitudes_matrix))) {
    rownames(amplitudes_matrix) <- design_info$conditions
  }
  
  # Create result object
  result <- structure(
    list(
      hrf_shapes = hrf_shapes,
      amplitudes = amplitudes_matrix,
      trial_amplitudes = trial_amplitudes,
      manifold_coords = smoothing_result$Xi_smooth,
      qc_metrics = qc_metrics,
      metadata = metadata,
      design_matrices = X_condition_list,
      voxel_indices = voxel_indices,
      data_info = data_info,
      call = mc
    ),
    class = c("mhrf_result", "list")
  )
  
  # Save if requested
  if (save_intermediate) {
    save_file <- file.path(output_dir, "mhrf_result.rds")
    saveRDS(result, save_file)
    progress$message(sprintf("Results saved to: %s", save_file))
  }
  
  # Generate report if requested
  if (params$generate_qc_report) {
    progress$update("Generating QC report...")
    report_file <- generate_qc_report(
      result,
      output_dir = output_dir,
      open_report = interactive()
    )
    progress$message(sprintf("QC report saved to: %s", report_file))
  }
  
  progress$complete(sprintf("Analysis completed in %.1f seconds", 
                           metadata$runtime_seconds))
  
  return(result)
}


# Helper functions (internal) -------------------------------------------------

#' Create progress tracker
#' @keywords internal
.create_progress_tracker <- function(verbose_level) {
  self <- list(
    level = verbose_level,
    start_time = NULL
  )
  
  self$start <- function(task) {
    self$start_time <<- Sys.time()
    if (self$level >= 1) {
      cat("\n", task, "\n", sep = "")
      cat(rep("=", nchar(task)), "\n", sep = "")
    }
  }
  
  self$update <- function(message, level = 1) {
    if (self$level >= level) {
      if (level == 1) {
        cat("\n", message, "\n", sep = "")
      } else {
        cat("  ", message, "\n", sep = "")
      }
    }
  }
  
  self$message <- function(message) {
    if (self$level >= 1) {
      cat("  ℹ ", message, "\n", sep = "")
    }
  }
  
  self$complete <- function(message = NULL) {
    if (self$level >= 1) {
      if (!is.null(message)) {
        cat("\n✓ ", message, "\n", sep = "")
      } else {
        cat("\n✓ Complete\n")
      }
    }
  }
  
  return(self)
}


#' Prepare data inputs
#' @keywords internal
.prepare_data_inputs <- function(Y_data, events, TR, voxel_mask, progress) {
  
  # Handle different input types
  if (inherits(Y_data, "matrix")) {
    Y_matrix <- Y_data
    voxel_coords <- NULL
    
  } else if (inherits(Y_data, c("NeuroVec", "NeuroVol"))) {
    # Handle neuroim2 objects
    progress$message("Converting neuroimaging data to matrix format")
    
    if (!requireNamespace("neuroim2", quietly = TRUE)) {
      stop("Package 'neuroim2' required for neuroimaging data. Please install it.")
    }
    
    # Extract data and coordinates
    Y_matrix <- as.matrix(Y_data)
    voxel_coords <- coordinates(Y_data)
    
  } else {
    stop("Y_data must be a matrix or NeuroVec/NeuroVol object")
  }
  
  # Get dimensions
  n_timepoints <- nrow(Y_matrix)
  n_voxels_total <- ncol(Y_matrix)
  
  # Apply mask if provided
  if (!is.null(voxel_mask)) {
    if (is.logical(voxel_mask)) {
      keep_voxels <- which(voxel_mask)
    } else if (is.numeric(voxel_mask)) {
      keep_voxels <- which(voxel_mask > 0)
    } else {
      stop("voxel_mask must be logical or numeric")
    }
    
    Y_matrix <- Y_matrix[, keep_voxels, drop = FALSE]
    if (!is.null(voxel_coords)) {
      voxel_coords <- voxel_coords[keep_voxels, , drop = FALSE]
    }
    voxel_indices <- keep_voxels
  } else {
    voxel_indices <- 1:n_voxels_total
  }
  
  n_voxels <- ncol(Y_matrix)
  
  progress$message(sprintf("Data dimensions: %d timepoints x %d voxels", 
                          n_timepoints, n_voxels))
  
  return(list(
    Y_matrix = Y_matrix,
    n_timepoints = n_timepoints,
    n_voxels = n_voxels,
    n_voxels_original = n_voxels_total,
    voxel_indices = voxel_indices,
    voxel_coords = voxel_coords
  ))
}


#' Determine number of jobs
#' @keywords internal
.determine_n_jobs <- function(n_jobs) {
  if (n_jobs == -1) {
    parallel::detectCores()
  } else {
    min(n_jobs, parallel::detectCores())
  }
}


#' Create design matrices from events
#' @keywords internal
.create_design_matrices <- function(events, n_timepoints, TR, hrf_length = 25) {

  if (!is.data.frame(events)) {
    stop("events must be a data frame")
  }

  required_cols <- c("onset", "condition")
  missing_cols <- setdiff(required_cols, names(events))
  if (length(missing_cols) > 0) {
    stop("events data frame missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  if (!"duration" %in% names(events)) {
    events$duration <- 0
  }

  events$condition <- as.factor(events$condition)
  conditions <- levels(events$condition)
  n_conditions <- length(conditions)
  n_trials <- nrow(events)

  # Use fmrireg to generate raw design matrices
  sframe <- fmrireg::sampling_frame(blocklens = n_timepoints, TR = TR)
  raw_basis <- HRF_RAW_EVENT_BASIS(hrf_length, TR)

  ev_model <- fmrireg::event_model(
    formula = ~ hrf(onset, basis = raw_basis, by = condition),
    data = events,
    sampling_frame = sframe,
    drop_empty = TRUE
  )

  X_full <- fmrireg::design_matrix(ev_model)

  term_tag <- names(ev_model$terms)[1]
  X_condition_list <- vector("list", n_conditions)
  for (i in seq_along(conditions)) {
    token <- fmrireg:::level_token("condition", conditions[i])
    prefix <- paste0(term_tag, "_", token)
    cols <- grep(paste0("^", prefix, "_b"), colnames(X_full))
    X_condition_list[[i]] <- X_full[, cols, drop = FALSE]
  }
  names(X_condition_list) <- conditions

  # Trial-wise design matrices using regressor evaluation
  times <- sframe$time
  X_trial_list <- vector("list", n_trials)
  for (j in seq_len(n_trials)) {
    reg <- fmrireg::regressor(
      onsets = events$onset[j],
      hrf = raw_basis,
      duration = events$duration[j],
      amplitude = 1,
      span = raw_basis$span
    )
    vals <- fmrireg::evaluate(reg, times)
    X_trial_list[[j]] <- matrix(vals, ncol = hrf_length)
  }

  return(list(
    X_condition_list = X_condition_list,
    X_trial_list = X_trial_list,
    n_conditions = n_conditions,
    n_trials = n_trials,
    conditions = conditions
  ))
}


#' Prepare HRF library
#' @keywords internal
.prepare_hrf_library <- function(hrf_library, p_hrf, TR, params) {
  time_points <- seq(0, by = TR, length.out = p_hrf)

  if (inherits(hrf_library, "HRF")) {
    # Single fmrireg HRF object
    L_library <- matrix(as.numeric(fmrireg::evaluate(hrf_library, time_points)),
                        ncol = 1)
    n_hrfs <- 1L
    library_type <- "hrf_object"

  } else if (is.list(hrf_library) && all(sapply(hrf_library, inherits, "HRF"))) {
    # List of fmrireg HRF objects
    L_library <- do.call(cbind, lapply(hrf_library, function(h) {
      as.numeric(fmrireg::evaluate(h, time_points))
    }))
    n_hrfs <- length(hrf_library)
    library_type <- "hrf_list"

  } else if (is.character(hrf_library)) {
    # Predefined library names
    if (hrf_library %in% c("auto", "gamma", "gamma_grid")) {
      hrf_objs <- create_gamma_grid_library(TR_precision = TR,
                                            hrf_duration = TR * (p_hrf - 1))
    } else if (hrf_library == "spmg1") {
      hrf_objs <- list(
        fmrireg::HRF_SPMG1,
        fmrireg::HRF_SPMG2,
        fmrireg::HRF_SPMG3
      )
    } else if (hrf_library == "flobs") {
      hrf_objs <- create_flobs_library(TR_precision = TR,
                                       hrf_duration = TR * (p_hrf - 1))
    } else {
      stop("Unknown HRF library type: ", hrf_library)
    }

    L_library <- do.call(cbind, lapply(hrf_objs, function(h) {
      as.numeric(fmrireg::evaluate(h, time_points))
    }))
    n_hrfs <- length(hrf_objs)
    library_type <- hrf_library

  } else if (is.matrix(hrf_library)) {
    # User-provided matrix
    L_library <- hrf_library
    n_hrfs <- ncol(L_library)
    library_type <- "custom"

  } else {
    stop("hrf_library must be a supported string, HRF object, list of HRF objects, or matrix")
  }
  
  # Quality check
  lib_quality <- check_hrf_library_quality(L_library)
  if (!lib_quality$is_good_quality) {
    warning("HRF library has quality issues. Consider using a different library.")
  }
  
  return(list(
    L_library = L_library,
    n_hrfs = n_hrfs,
    library_type = library_type,
    quality = lib_quality
  ))
}


#' Run manifold construction
#' @keywords internal
.run_manifold_construction <- function(L_library, params, progress) {
  
  # Calculate affinity matrix
  S_markov <- calculate_manifold_affinity_core(
    L_library_matrix = L_library,
    k_local_nn_for_sigma = params$k_local_nn_for_sigma
  )
  
  # Get manifold basis
  if (params$use_robust_manifold) {
    manifold <- get_manifold_basis_reconstructor_robust(
      S_markov_matrix = S_markov,
      L_library_matrix = L_library,
      m_manifold_dim_target = params$m_manifold_dim_target,
      m_manifold_dim_min_variance = params$m_manifold_dim_min_variance,
      fallback_to_pca = params$fallback_to_pca
    )
  } else {
    manifold <- get_manifold_basis_reconstructor_core(
      S_markov_matrix = S_markov,
      L_library_matrix = L_library,
      m_manifold_dim_target = params$m_manifold_dim_target,
      m_manifold_dim_min_variance = params$m_manifold_dim_min_variance
    )
  }
  
  progress$message(sprintf("Manifold constructed: %d dimensions (method: %s)",
                          manifold$m_final_dim, 
                          manifold$method_used %||% "diffusion_map"))
  
  return(list(
    B_reconstructor = manifold$B_reconstructor_matrix,
    Phi_coords = manifold$Phi_coords_matrix,
    eigenvalues = manifold$eigenvalues_S_vector,
    m_final = manifold$m_final_dim,
    method_used = manifold$method_used %||% "diffusion_map"
  ))
}


#' Run voxelwise estimation
#' @keywords internal  
.run_voxelwise_estimation <- function(Y_data, X_condition_list, manifold, 
                                     params, progress) {
  
  # Transform designs to manifold basis
  XB_list <- transform_designs_to_manifold_basis_core(
    X_condition_list_proj_matrices = X_condition_list,
    B_reconstructor_matrix = manifold$B_reconstructor
  )
  
  # Solve for gamma coefficients
  Gamma <- solve_glm_for_gamma_core(
    Z_list_of_matrices = XB_list,
    Y_proj_matrix = Y_data,
    lambda_gamma = params$lambda_gamma,
    orthogonal_approx_flag = params$orthogonal_approx %||% FALSE
  )
  
  # Extract Xi and Beta
  if (params$use_robust_svd) {
    svd_result <- extract_xi_beta_raw_svd_robust(
      Gamma_coeffs_matrix = Gamma,
      m_manifold_dim = manifold$m_final,
      k_conditions = length(X_condition_list)
    )
  } else {
    svd_result <- extract_xi_beta_raw_svd_core(
      Gamma_coeffs_matrix = Gamma,
      m_manifold_dim = manifold$m_final,
      k_conditions = length(X_condition_list)
    )
  }
  
  # Apply identifiability constraints
  # Note: The current implementation requires additional parameters
  # For now, we'll use the raw results
  Xi_ident <- svd_result$Xi_raw_matrix
  Beta_ident <- svd_result$Beta_raw_matrix
  
  # Basic sign correction - ensure positive mean amplitude
  for (v in 1:ncol(Beta_ident)) {
    if (mean(Beta_ident[, v]) < 0) {
      Beta_ident[, v] <- -Beta_ident[, v]
      Xi_ident[, v] <- -Xi_ident[, v]
    }
  }
  
  return(list(
    Xi_ident = Xi_ident,
    Beta_ident = Beta_ident,
    Gamma = Gamma
  ))
}


#' Run spatial smoothing
#'
#' @param Xi_matrix m x V matrix of manifold coordinates
#' @param voxel_coords V x 3 matrix of voxel coordinates
#' @param params List of smoothing parameters
#' @param progress Progress reporter object
#' @param Y_data Optional n x V data matrix for computing local SNR
#' @keywords internal
.run_spatial_smoothing <- function(Xi_matrix, voxel_coords, params, progress,
                                   Y_data = NULL) {
  
  n_voxels <- ncol(Xi_matrix)
  
  # Create spatial graph if coordinates available
  if (!is.null(voxel_coords) && nrow(voxel_coords) == n_voxels) {
    # Ensure 3D coordinates
    if (ncol(voxel_coords) == 2) {
      voxel_coords <- cbind(voxel_coords, rep(1, n_voxels))
    }
    
    # Create graph Laplacian
    L_spatial <- make_voxel_graph_laplacian_core(
      voxel_coords = voxel_coords,
      num_neighbors = params$num_neighbors_Lsp
    )
    
  } else {
    # No spatial information - use identity (no smoothing)
    progress$message("No spatial coordinates available - skipping spatial smoothing")
    L_spatial <- Matrix::Diagonal(n_voxels)
  }
  
  # Apply smoothing
  if (params$adaptive_smoothing && exists("compute_local_snr")) {
    # Adaptive smoothing based on SNR
    local_snr <- if (!is.null(Y_data)) {
      compute_local_snr(Y_data)
    } else {
      rep(1, n_voxels)
    }

    Xi_smooth <- apply_spatial_smoothing_adaptive(
      Xi_ident_matrix = Xi_matrix,
      L_sp_sparse_matrix = L_spatial,
      lambda_spatial_smooth = params$lambda_spatial_smooth,
      local_snr = local_snr,
      edge_preserve = params$edge_preserve %||% FALSE,
      voxel_coords = voxel_coords
    )
  } else {
    # Standard smoothing
    Xi_smooth <- apply_spatial_smoothing_core(
      Xi_ident_matrix = Xi_matrix,
      L_sp_sparse_matrix = L_spatial,
      lambda_spatial_smooth = params$lambda_spatial_smooth
    )
  }
  
  return(list(
    Xi_smooth = Xi_smooth,
    L_spatial = L_spatial
  ))
}


#' Run trial estimation
#' @keywords internal
.run_trial_estimation <- function(Y_data, X_trial_list, hrf_shapes, 
                                 params, progress) {
  
  n_voxels <- ncol(Y_data)
  n_trials <- length(X_trial_list)
  
  # Initialize storage
  trial_amplitudes <- matrix(NA, n_trials, n_voxels)
  
  # Process in chunks for memory efficiency
  chunk_size <- params$chunk_size %||% 1000
  n_chunks <- ceiling(n_voxels / chunk_size)
  
  if (params$n_jobs > 1 && requireNamespace("future", quietly = TRUE)) {
    # Parallel processing
    progress$message(sprintf("Processing %d voxels in parallel (%d jobs)", 
                            n_voxels, params$n_jobs))
    
    # Set up parallel backend
    oplan <- future::plan(future::multisession, workers = params$n_jobs)
    on.exit(future::plan(oplan))
    
    # Process chunks in parallel
    chunk_results <- future.apply::future_lapply(1:n_chunks, function(chunk) {
      start_idx <- (chunk - 1) * chunk_size + 1
      end_idx <- min(chunk * chunk_size, n_voxels)
      chunk_voxels <- start_idx:end_idx
      
      chunk_result <- matrix(NA, n_trials, length(chunk_voxels))
      
      for (v_idx in seq_along(chunk_voxels)) {
        v <- chunk_voxels[v_idx]
        
        lss_result <- run_lss_for_voxel_corrected(
          y_voxel = Y_data[, v],
          X_trial_list = X_trial_list,
          h_voxel = hrf_shapes[, v],
          TR = params$TR,
          lambda = params$lambda_ridge_Alss
        )
        
        chunk_result[, v_idx] <- lss_result$beta_trials
      }
      
      return(chunk_result)
    })
    
    # Combine results
    for (chunk in 1:n_chunks) {
      start_idx <- (chunk - 1) * chunk_size + 1
      end_idx <- min(chunk * chunk_size, n_voxels)
      trial_amplitudes[, start_idx:end_idx] <- chunk_results[[chunk]]
    }
    
  } else {
    # Sequential processing
    if (params$verbose_level >= 2) {
      pb <- create_progress_bar(n_voxels)
    }
    
    for (v in 1:n_voxels) {
      lss_result <- run_lss_for_voxel_corrected(
        y_voxel = Y_data[, v],
        X_trial_list = X_trial_list,
        h_voxel = hrf_shapes[, v],
        TR = params$TR,
        lambda = params$lambda_ridge_Alss
      )
      
      trial_amplitudes[, v] <- lss_result$beta_trials
      
      if (params$verbose_level >= 2 && v %% 100 == 0) {
        pb <- update_progress_bar(pb, 100)
      }
    }
  }
  
  return(trial_amplitudes)
}


#' Compute QC metrics
#'
#' Estimates basic quality metrics for the fitted model.  A sample of voxels is
#' reconstructed using the estimated HRFs and amplitudes and the resulting
#' R\eqn{^2} statistics are summarized.
#'
#' @param Y_data n \times V matrix of observed time series
#' @param X_condition_list list of condition design matrices
#' @param hrf_shapes p \times V matrix of estimated HRF shapes
#' @param amplitudes k \times V matrix of condition amplitudes
#' @param manifold_coords m \times V matrix of manifold coordinates
#' @param params List of pipeline parameters
#' @keywords internal
.compute_qc_metrics <- function(Y_data, X_condition_list, hrf_shapes, amplitudes,
                               manifold_coords, params) {
  
  # Basic metrics
  qc <- list(
    n_voxels_analyzed = ncol(Y_data),
    mean_amplitude = mean(amplitudes),
    sd_amplitude = sd(amplitudes),
    percent_negative_amp = 100 * mean(amplitudes < 0),
    manifold_variance = apply(manifold_coords, 1, var)
  )
  
  # HRF statistics
  hrf_stats <- extract_hrf_stats(hrf_shapes, TR = params$TR)
  qc$hrf_stats <- hrf_stats
  
  # Reconstruction quality (sample)
  if (ncol(Y_data) > 100) {
    # Sample voxels for efficiency
    sample_voxels <- sample(ncol(Y_data), 100)
  } else {
    sample_voxels <- 1:ncol(Y_data)
  }
  
  # Reconstruct a subset of voxels to evaluate fit quality
  r2_voxels <- numeric(length(sample_voxels))
  for (i in seq_along(sample_voxels)) {
    v <- sample_voxels[i]
    y <- Y_data[, v]

    # Predicted time series from condition regressors
    y_pred <- rep(0, nrow(Y_data))
    h_v <- hrf_shapes[, v]
    for (j in seq_along(X_condition_list)) {
      y_pred <- y_pred + (X_condition_list[[j]] %*% h_v) * amplitudes[j, v]
    }

    resid <- y - y_pred
    denom <- sum((y - mean(y))^2)
    if (denom > 0) {
      r2_voxels[i] <- 1 - sum(resid^2) / denom
    } else {
      r2_voxels[i] <- NA_real_
    }
  }

  qc$mean_r_squared <- mean(r2_voxels, na.rm = TRUE)
  qc$r_squared_voxels <- r2_voxels
  qc$diagnostics <- list(r2_voxelwise = r2_voxels)
  qc$quality_flags <- create_qc_flags(qc)
  
  return(qc)
}


# Additional helper utilities -------------------------------------------------

#' Get coordinate extraction function based on data type
#' @keywords internal
coordinates <- function(x) {
  UseMethod("coordinates")
}

#' Default coordinates method
#' @keywords internal
coordinates.default <- function(x) {
  NULL
}

#' Extract voxel coordinates from NeuroVol
#' @keywords internal  
coordinates.NeuroVol <- function(x) {
  if (requireNamespace("neuroim2", quietly = TRUE)) {
    # Extract 3D coordinates
    coords <- neuroim2::coords(x)
    return(as.matrix(coords))
  }
  return(NULL)
}

#' Extract voxel coordinates from NeuroVec
#' @keywords internal
coordinates.NeuroVec <- function(x) {
  if (requireNamespace("neuroim2", quietly = TRUE)) {
    # Extract coordinates from vector
    space_info <- neuroim2::space(x)
    indices <- neuroim2::indices(x)
    
    # Convert indices to coordinates
    coords <- neuroim2::index_to_coord(space_info, indices)
    return(as.matrix(coords))
  }
  return(NULL)
}

#' Create progress bar for sequential processing
#' @keywords internal
create_progress_bar <- function(total) {
  list(
    total = total,
    current = 0,
    start_time = Sys.time(),
    pb = txtProgressBar(min = 0, max = total, style = 3)
  )
}

#' Update progress bar
#' @keywords internal
update_progress_bar <- function(pb_obj, increment = 1) {
  pb_obj$current <- pb_obj$current + increment
  setTxtProgressBar(pb_obj$pb, pb_obj$current)
  
  # Estimate time remaining
  elapsed <- as.numeric(difftime(Sys.time(), pb_obj$start_time, units = "secs"))
  if (pb_obj$current > 0) {
    rate <- pb_obj$current / elapsed
    remaining <- (pb_obj$total - pb_obj$current) / rate
    if (remaining > 60) {
      message(sprintf("\nEstimated time remaining: %.1f minutes", remaining / 60))
    }
  }
  
  # Close if complete
  if (pb_obj$current >= pb_obj$total) {
    close(pb_obj$pb)
  }
  
  return(pb_obj)
}
