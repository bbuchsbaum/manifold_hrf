# User-Facing Interface for M-HRF-LSS
# Integrates with fmrireg infrastructure

#' Fit M-HRF-LSS Model
#'
#' Main interface for manifold-guided HRF estimation and trial-wise deconvolution.
#' Follows fmrireg patterns but extends with manifold-based HRF estimation.
#'
#' @param formula Formula specification for event model (same as fmrireg)
#' @param dataset An fmri_dataset object from fmrireg
#' @param baseline_model Optional baseline model from fmrireg::baseline_model()
#' @param hrf_library Source for HRF library. Can be:
#'   - Character: "canonical" (uses fmrireg HRFs), "flobs", "gamma_grid"
#'   - List of fmrireg HRF objects (e.g., list(HRF_SPMG1, HRF_GAMMA))
#'   - Matrix of HRF shapes (p x N)
#'   - Path to saved HRF library
#' @param manifold_params List of manifold construction parameters or preset name
#' @param estimation What to estimate: "condition", "trial", or "both"
#' @param spatial_info Optional NeuroSpace or voxel coordinates for spatial smoothing
#' @param strategy Processing strategy: "global", "runwise", "chunked"
#' @param nchunks Number of chunks for memory efficiency
#' @param robust Use robust estimation (downweight outliers)
#' @param progress Show progress bar
#' @param verbose Print detailed messages
#' @param ... Additional arguments
#'
#' @return An \code{mhrf_lss_result} object with components:
#'   \describe{
#'     \item{parameters}{List of parameter settings used}
#'     \item{manifold}{HRF manifold object}
#'     \item{hrf}{List of raw and smoothed HRF matrices}
#'     \item{beta}{List of condition and trial beta estimates}
#'     \item{qc}{Quality control metrics}
#'     \item{diagnostics}{Additional diagnostic information}
#'   }
#'
#' @examples
#' \dontrun{
#' # Load fmrireg for dataset creation
#' library(fmrireg)
#' 
#' # Create dataset (using fmrireg)
#' dset <- fmri_dataset(
#'   scans = "bold.nii.gz",
#'   mask = "mask.nii.gz", 
#'   TR = 2,
#'   run_length = c(200, 200)
#' )
#' 
#' # Basic M-HRF-LSS fit
#' fit <- mhrf_lss(
#'   ~ hrf(condition),
#'   dataset = dset,
#'   hrf_library = "canonical",
#'   estimation = "both"
#' )
#' 
#' # Advanced with custom parameters
#' fit_custom <- mhrf_lss(
#'   ~ hrf(cond1) + hrf(cond2),
#'   dataset = dset,
#'   baseline_model = baseline_model(degree = 2),
#'   hrf_library = list(HRF_SPMG1, HRF_SPMG2, HRF_SPMG3, HRF_GAMMA),
#'   manifold_params = "balanced",  # or custom list
#'   spatial_info = dset$mask,
#'   robust = TRUE
#' )
#'
#' # Access estimated HRFs
#' head(fit$hrf$smoothed)
#' }
#'
#' @description Deprecated wrapper for `mhrf_analyze`. Use that function
#'   in new code.
#' @keywords internal
mhrf_lss <- function(formula,
                     dataset,
                     baseline_model = NULL,
                     hrf_library = "canonical",
                     manifold_params = "balanced",
                     estimation = c("condition", "trial", "both"),
                     spatial_info = NULL,
                     strategy = c("global", "runwise", "chunked"),
                     nchunks = 10,
                     robust = FALSE,
                     progress = TRUE,
                     verbose = TRUE,
                     ...) {

  .Deprecated("mhrf_analyze")
  

  
  estimation <- match.arg(estimation)
  strategy <- match.arg(strategy)
  
  # Extract data dimensions
  data_mat <- fmrireg::get_data_matrix(dataset)
  n_time <- nrow(data_mat)
  n_voxels <- ncol(data_mat)

  # Basic validation using helper functions
  .validate_Y_data(data_mat)
  .validate_events(dataset$event_table, n_time, dataset$TR)
  
  if (verbose) {
    message(sprintf("M-HRF-LSS: %d timepoints, %d voxels", n_time, n_voxels))
  }
  
  # Step 1: Build event model using fmrireg
  if (verbose) message("Building event model...")
  
  # Get sampling frame from dataset
  sframe <- fmrireg::sampling_frame(
    blocklens = dataset$run_length,
    TR = dataset$TR
  )
  
  # Create event model
  event_mod <- fmrireg::event_model(
    formula = formula,
    data = dataset$event_table,
    block = ~ run,  # Assuming 'run' column exists
    sampling_frame = sframe,
    drop_empty = TRUE
  )
  
  # Extract design matrices for conditions and trials
  design_info <- extract_design_info(event_mod, sframe)
  validate_design_matrix_list(design_info$X_condition_list, n_time)
  
  # Step 2: Create HRF manifold
  if (verbose) message("Constructing HRF manifold...")
  
  manifold <- create_hrf_manifold(
    hrf_library = hrf_library,
    params = manifold_params,
    TR = dataset$TR,
    verbose = verbose
  )
  
  # Step 3: Prepare spatial information if provided
  voxel_coords <- NULL
  if (!is.null(spatial_info)) {
    voxel_coords <- extract_voxel_coordinates(spatial_info, dataset$mask)
  }
  
  # Step 4: Setup baseline/confounds
  if (is.null(baseline_model)) {
    # Default baseline: intercept per run + linear drift
    baseline_model <- fmrireg::baseline_model(
      basis = "poly",
      degree = 1,
      sframe = sframe
    )
  }
  
  # Get confound matrix
  Z_confounds <- fmrireg::design_matrix(baseline_model)
  validate_confounds_matrix(Z_confounds, n_time)
  
  # Step 5: Check for outliers if robust
  outlier_weights <- NULL
  if (robust) {
    if (verbose) message("Detecting outliers...")
    outlier_weights <- detect_outlier_timepoints(data_mat)
  }
  
  # Step 6: Run M-HRF-LSS pipeline
  if (verbose) message("Running M-HRF-LSS estimation...")
  
  # Determine processing strategy
  if (strategy == "chunked" || n_voxels > 50000) {
    # Use chunked processing for large datasets
    results <- run_mhrf_lss_chunked(
      Y_data = data_mat,
      design_info = design_info,
      manifold = manifold,
      Z_confounds = Z_confounds,
      voxel_coords = voxel_coords,
      params = manifold$parameters,
      outlier_weights = outlier_weights,
      estimation = estimation,
      nchunks = nchunks,
      progress = progress
    )
  } else {
    # Standard processing
    results <- run_mhrf_lss_standard(
      Y_data = data_mat,
      design_info = design_info,
      manifold = manifold,
      Z_confounds = Z_confounds,
      voxel_coords = voxel_coords,
      params = manifold$parameters,
      outlier_weights = outlier_weights,
      estimation = estimation,
      progress = progress
    )
  }
  
  # Step 7: Package results
  if (verbose) message("Packaging results...")
  
  fit <- structure(
    list(
      parameters = list(
        manifold = manifold$parameters,
        robust = robust
      ),
      manifold = manifold,
      hrf = list(
        raw = results$H_raw,
        smoothed = results$H_shapes
      ),
      beta = list(
        condition_initial = results$Beta_condition_initial,
        condition_final = results$Beta_condition,
        trial = if (estimation %in% c("trial", "both")) results$Beta_trial else NULL
      ),
      qc = results$diagnostics,
      diagnostics = results$diagnostics,
      call = match.call()
    ),
    class = "mhrf_lss_result"
  )
  
  if (verbose) message("M-HRF-LSS fitting complete!")
  
  return(fit)
}


#' Create HRF Manifold from Various Sources
#'
#' Internal function to create HRF manifold compatible with fmrireg
#'
#' @keywords internal
create_hrf_manifold <- function(hrf_library, params, TR, verbose = TRUE) {
  
  # Handle parameter presets
  if (is.character(params) && length(params) == 1) {
    if (verbose) {
      params <- get_preset_params(params)
    } else {
      params <- suppressMessages(get_preset_params(params))
    }
  }
  
  # Determine sampling grid for evaluating HRFs
  hrf_duration <- params$hrf_duration %||% 24
  time_points <- seq(0, hrf_duration, by = TR)

  if (inherits(hrf_library, "HRF")) {
    L_library <- matrix(as.numeric(fmrireg::evaluate(hrf_library, time_points)),
                        ncol = 1)
  } else if (is.list(hrf_library) && all(sapply(hrf_library, inherits, "HRF"))) {
    L_library <- do.call(cbind, lapply(hrf_library, function(h) {
      as.numeric(fmrireg::evaluate(h, time_points))
    }))
  } else if (is.character(hrf_library) && length(hrf_library) == 1) {
    if (hrf_library %in% c("canonical", "spmg1")) {
      hrf_objs <- list(
        fmrireg::HRF_SPMG1,
        fmrireg::HRF_SPMG2,
        fmrireg::HRF_SPMG3
      )
    } else if (hrf_library %in% c("gamma", "gamma_grid", "auto")) {
      hrf_objs <- create_gamma_grid_library(TR_precision = TR,
                                            hrf_duration = hrf_duration)
    } else if (hrf_library == "flobs") {
      hrf_objs <- create_flobs_library(TR_precision = TR,
                                       hrf_duration = hrf_duration)
    } else {
      stop("Unknown HRF library type: ", hrf_library)
    }

    L_library <- do.call(cbind, lapply(hrf_objs, function(h) {
      as.numeric(fmrireg::evaluate(h, time_points))
    }))
  } else if (is.matrix(hrf_library)) {
    L_library <- hrf_library
  } else {
    stop("Invalid hrf_library format")
  }
  
  # Calculate affinity matrix
  S_markov <- calculate_manifold_affinity_core(
    L_library_matrix = L_library,
    k_local_nn_for_sigma = params$k_local_nn_for_sigma %||% 10
  )
  
  # Create manifold with robust construction
  if (verbose) {
    manifold <- get_manifold_basis_reconstructor_robust(
      S_markov_matrix = S_markov,
      L_library_matrix = L_library,
      m_manifold_dim_target = params$m_manifold_dim_target,
      m_manifold_dim_min_variance = params$m_manifold_dim_min_variance %||% 0.95,
      fallback_to_pca = TRUE
    )
  } else {
    manifold <- withCallingHandlers({
      get_manifold_basis_reconstructor_robust(
        S_markov_matrix = S_markov,
        L_library_matrix = L_library,
        m_manifold_dim_target = params$m_manifold_dim_target,
        m_manifold_dim_min_variance = params$m_manifold_dim_min_variance %||% 0.95,
        fallback_to_pca = TRUE
      )
    }, warning = function(w) {
      if (verbose) message("Warning in manifold construction: ", conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  }
  
  # Add additional parameters for compatibility
  manifold$parameters <- params
  manifold$library_hrfs <- L_library
  manifold$m_manifold_dim <- manifold$m_final_dim %||% manifold$m_manifold_dim
  
  return(manifold)
}


#' Extract Design Information from fmrireg Event Model
#'
#' @keywords internal
extract_design_info <- function(event_model, sframe) {
  
  # Get condition-wise design matrices
  term_mats <- fmrireg::term_matrices(event_model)
  
  # Extract individual trial information
  event_tab <- fmrireg::event_table(event_model)
  
  # Create condition design matrices (for standard GLM)
  # Convert term_mats to list of matrices if needed
  if (is.matrix(term_mats)) {
    # If term_mats is a single matrix, wrap in a list
    X_condition_list <- list(term_mats)
  } else if (is.list(term_mats)) {
    # If it's already a list, ensure each element is a matrix
    X_condition_list <- lapply(term_mats, function(x) {
      if (is.matrix(x)) {
        return(x)
      } else {
        # Try to convert to matrix
        return(as.matrix(x))
      }
    })
  } else {
    # Fallback: try to convert to matrix and wrap in list
    X_condition_list <- list(as.matrix(term_mats))
  }
  
  # Create trial-wise design matrices (for LSS)
  # This requires reconstructing individual trial regressors
  X_trial_list <- create_trial_matrices(event_model, event_tab, sframe)
  
  return(list(
    X_condition_list = X_condition_list,
    X_trial_list = X_trial_list,
    event_table = event_tab,
    n_conditions = length(X_condition_list),
    n_trials = nrow(event_tab)
  ))
}


#' Create Trial-wise Design Matrices
#'
#' @keywords internal
create_trial_matrices <- function(event_model, event_table, sframe) {
  
  n_trials <- nrow(event_table)
  n_time <- nrow(sframe)
  
  # Extract HRF from the event model's first term
  # In fmrireg, HRF is usually stored in the term specification
  terms <- fmrireg::terms(event_model)
  
  # Get the HRF object from the first term (assuming all use same HRF)
  if (length(terms) > 0 && "hrf" %in% names(terms[[1]])) {
    hrf_obj <- terms[[1]]$hrf
  } else {
    # Default to canonical HRF
    hrf_obj <- fmrireg::HRF_SPMG1
  }
  
  p <- hrf_obj$nbasis
  
  # Get the sampling times
  times <- sframe$time
  
  X_trial_list <- list()
  
  for (i in 1:n_trials) {
    # Create single trial design matrix
    onset_time <- event_table$onset[i]
    duration <- event_table$duration[i] %||% 0
    block_id <- event_table$block[i] %||% 1
    
    # Find which block this trial belongs to
    block_start_idx <- which(sframe$block == block_id)[1]
    block_times <- times[sframe$block == block_id]
    
    # Create regressor for this trial
    trial_reg <- fmrireg::regressor(
      onsets = onset_time,
      hrf = hrf_obj,
      duration = duration,
      amplitude = 1,
      span = hrf_obj$span %||% 24
    )
    
    # Evaluate at the block times
    hrf_values <- fmrireg::evaluate(trial_reg, block_times - block_times[1])
    
    # Create full design matrix (accounting for all timepoints)
    X_trial <- matrix(0, n_time, p)
    block_indices <- which(sframe$block == block_id)
    
    if (length(hrf_values) == length(block_indices) * p) {
      # Reshape if multiple basis functions
      X_trial[block_indices, ] <- matrix(hrf_values, ncol = p)
    } else if (length(hrf_values) == length(block_indices)) {
      # Single basis function
      X_trial[block_indices, 1] <- hrf_values
    }
    
    X_trial_list[[i]] <- X_trial
  }
  
  return(X_trial_list)
}


#' Run Standard M-HRF-LSS Pipeline
#'
#' @keywords internal
run_mhrf_lss_standard <- function(Y_data, design_info, manifold, Z_confounds,
                                  voxel_coords, params, outlier_weights,
                                  estimation, progress) {
  
  # This wraps our core functions in the right sequence
  
  # 1. Project out confounds
  proj_result <- project_out_confounds_core(
    Y_data_matrix = Y_data,
    X_list_of_matrices = design_info$X_condition_list,
    Z_confounds_matrix = Z_confounds
  )
  
  # 2. Transform to manifold basis
  Z_list <- transform_designs_to_manifold_basis_core(
    X_condition_list_proj_matrices = proj_result$X_list_proj_matrices,
    B_reconstructor_matrix = manifold$B_reconstructor_matrix
  )
  
  # 3. Solve for gamma coefficients
  # If using outlier weights, modify Y_proj
  Y_proj_weighted <- proj_result$Y_proj_matrix
  if (!is.null(outlier_weights)) {
    Y_proj_weighted <- Y_proj_weighted * sqrt(outlier_weights)
  }
  
  Gamma_coeffs <- solve_glm_for_gamma_core(
    Z_list_of_matrices = Z_list,
    Y_proj_matrix = Y_proj_weighted,
    lambda_gamma = params$lambda_gamma
  )
  
  # 4. Extract Xi and Beta
  xi_beta <- extract_xi_beta_raw_svd_robust(
    Gamma_coeffs_matrix = Gamma_coeffs,
    m_manifold_dim = manifold$m_manifold_dim,
    k_conditions = design_info$n_conditions
  )
  
  # 5. Apply identifiability constraints
  # Use canonical HRF as reference
  h_ref <- manifold$library_hrfs[, 1]  # First HRF as reference
  
  # Ensure h_ref is a numeric vector
  if (is.matrix(h_ref)) {
    h_ref <- as.numeric(h_ref)
  }
  
  ident_result <- apply_intrinsic_identifiability_core(
    Xi_raw_matrix = xi_beta$Xi_raw_matrix,
    Beta_raw_matrix = xi_beta$Beta_raw_matrix,
    B_reconstructor_matrix = manifold$B_reconstructor_matrix,
    h_ref_shape_vector = h_ref,
    ident_scale_method = "l2_norm",
    ident_sign_method = params$ident_sign_method %||% "canonical_correlation",
    Y_proj_matrix = proj_result$Y_proj_matrix,
    X_condition_list_proj_matrices = proj_result$X_list_proj_matrices
  )

  # HRF shapes prior to spatial smoothing
  H_raw <- reconstruct_hrf_shapes_core(
    B_reconstructor_matrix = manifold$B_reconstructor_matrix,
    Xi_smoothed_matrix = ident_result$Xi_ident_matrix
  )

  # 6. Spatial smoothing if coordinates provided
  if (!is.null(voxel_coords)) {
    # Create spatial graph
    L_spatial <- make_voxel_graph_laplacian_core(
      voxel_coords_matrix = voxel_coords,
      num_neighbors_Lsp = params$num_neighbors_Lsp %||% 6
    )
    
    # Compute local SNR for adaptive smoothing
    local_snr <- compute_local_snr(Y_data, method = "temporal_variance")
    
    # Apply adaptive smoothing
    Xi_smoothed <- apply_spatial_smoothing_adaptive(
      Xi_ident_matrix = ident_result$Xi_ident_matrix,
      L_sp_sparse_matrix = L_spatial,
      lambda_spatial_smooth = params$lambda_spatial_smooth,
      local_snr = local_snr
    )
  } else {
    Xi_smoothed <- ident_result$Xi_ident_matrix
  }
  
  # 7. Reconstruct HRF shapes
  H_shapes <- reconstruct_hrf_shapes_core(
    B_reconstructor_matrix = manifold$B_reconstructor_matrix,
    Xi_smoothed_matrix = Xi_smoothed
  )
  
  # 8. Trial-wise estimation if requested
  Beta_trial <- NULL
  if (estimation %in% c("trial", "both")) {
    # Prepare LSS components
    lss_prep <- prepare_lss_fixed_components_core(
      A_lss_fixed_matrix = Z_confounds,
      intercept_col_index_in_Alss = 1,
      lambda_ridge_Alss = params$lambda_ridge_Alss %||% 1e-6
    )
    
    # Run LSS
    Beta_trial <- run_lss_voxel_loop_core(
      Y_proj_matrix = proj_result$Y_proj_matrix,
      X_trial_onset_list_of_matrices = design_info$X_trial_list,
      H_shapes_allvox_matrix = H_shapes,
      A_lss_fixed_matrix = Z_confounds,
      P_lss_matrix = lss_prep$P_lss_matrix,
      p_lss_vector = lss_prep$p_lss_vector,
      n_jobs = params$n_jobs %||% 1
    )
  }
  
  # 9. Re-estimate condition betas with final HRFs
  Beta_condition_final <- estimate_final_condition_betas_core(
    Y_proj_matrix = proj_result$Y_proj_matrix,
    X_condition_list_proj_matrices = proj_result$X_list_proj_matrices,
    H_shapes_allvox_matrix = H_shapes,
    lambda_beta_final = params$lambda_beta_final,
    control_alt_list = list(max_iter = 1)
  )
  
  # 10. Compute diagnostics
  diagnostics <- compute_pipeline_diagnostics(
    Y_data = Y_data,
    Y_proj = proj_result$Y_proj_matrix,
    H_shapes = H_shapes,
    Beta_condition = Beta_condition_final,
    Xi_smoothed = Xi_smoothed,
    TR_precision = params$TR_precision %||% 1
  )
  
  return(list(
    H_raw = H_raw,
    H_shapes = H_shapes,
    Xi_smoothed = Xi_smoothed,
    Beta_condition_initial = ident_result$Beta_ident_matrix,
    Beta_condition = Beta_condition_final,
    Beta_trial = Beta_trial,
    diagnostics = diagnostics
  ))
}


#' Compute Pipeline Diagnostics
#'
#' @param Y_data Original data matrix
#' @param Y_proj Confound-projected data matrix
#' @param H_shapes Matrix of reconstructed HRFs
#' @param Beta_condition Condition-level beta estimates
#' @param Xi_smoothed Smoothed manifold coordinates
#' @param TR_precision Sampling rate of the HRFs
#' @keywords internal
compute_pipeline_diagnostics <- function(Y_data, Y_proj, H_shapes,
                                       Beta_condition, Xi_smoothed,
                                       TR_precision = 1) {
  
  # This would compute R², HRF statistics, etc.
  # Simplified version:
  
  diagnostics <- list(
    n_voxels_processed = ncol(H_shapes),
    manifold_variance = apply(Xi_smoothed, 1, var),
    hrf_stats = extract_hrf_stats(H_shapes, TR_precision = TR_precision),
    timestamp = Sys.time()
  )
  
  return(diagnostics)
}


#' Extract Voxel Coordinates
#'
#' @param spatial_info Either a matrix of coordinates or a mask array
#' @param mask Optional mask to apply
#' @return Matrix of voxel coordinates or NULL
#' @keywords internal
extract_voxel_coordinates <- function(spatial_info, mask = NULL) {
  if (is.null(spatial_info)) {
    return(NULL)
  }
  
  if (is.matrix(spatial_info)) {
    # Already a coordinate matrix
    return(spatial_info)
  }
  
  if (is.array(spatial_info)) {
    # Extract coordinates from mask array
    if (!is.null(mask)) {
      # Use mask to select voxels
      indices <- which(mask != 0, arr.ind = TRUE)
      return(indices)
    } else {
      # Use all non-zero voxels
      indices <- which(spatial_info != 0, arr.ind = TRUE)
      return(indices)
    }
  }
  
  stop("spatial_info must be either a matrix or an array")
}


# S3 Methods for mhrf_lss_result

#' @export
print.mhrf_lss_result <- function(x, ...) {
  cat("M-HRF-LSS Result\n")
  cat("=================\n\n")

  if (!is.null(x$manifold)) {
    cat("Manifold dimensions:", x$manifold$m_manifold_dim, "\n")
  }

  if (!is.null(x$beta$condition_final)) {
    cat("Conditions:", nrow(x$beta$condition_final), "\n")
  }

  invisible(x)
}


#' Extract Coefficients from M-HRF-LSS Fit
#'
#' @export
coef.mhrf_lss_result <- function(object, type = c("condition", "trial", "hrf"), ...) {
  
  type <- match.arg(type)
  
  switch(type,
    condition = object$beta$condition_final,
    trial = object$beta$trial,
    hrf = object$hrf$smoothed
  )
}


#' Plot M-HRF-LSS Results
#'
#' @export
plot.mhrf_lss_result <- function(x, type = c("hrfs", "manifold", "diagnostics"),
                                 voxels = NULL, ...) {

  type <- match.arg(type)

  # Build a minimal object compatible with plotting helpers from
  # mhrf_result_methods.R
  plot_obj <- list(
    hrf_shapes = x$hrf$smoothed,
    amplitudes = x$beta$condition_final,
    manifold_coords = NULL,
    qc_metrics = x$qc,
    metadata = list(parameters = list(
      TR = x$parameters$manifold$TR_precision %||% 1
    ))
  )

  if (!is.null(x$diagnostics$Xi_smoothed)) {
    plot_obj$manifold_coords <- x$diagnostics$Xi_smoothed
  }

  if (type == "hrfs") {
    .plot_hrf_shapes(plot_obj, voxels = voxels, ...)

  } else if (type == "manifold") {
    if (is.null(plot_obj$manifold_coords)) {
      message("Manifold coordinates not available for plotting")
    } else {
      .plot_manifold_coords(plot_obj, ...)
    }

  } else if (type == "diagnostics") {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    .plot_hrf_shapes(plot_obj, voxels = voxels)
    .plot_amplitude_distribution(plot_obj)
    if (!is.null(plot_obj$manifold_coords)) {
      .plot_manifold_coords(plot_obj)
    } else {
      plot.new(); title("Manifold Coordinates")
      text(0.5, 0.5, "Not available")
    }
    .plot_quality_metrics(plot_obj)
  }

  invisible(x)
}


#' Run M-HRF-LSS Pipeline in Chunks
#'
#' Memory-efficient chunked processing for large datasets
#'
#' @keywords internal
run_mhrf_lss_chunked <- function(Y_data, design_info, manifold, Z_confounds,
                                voxel_coords, params, outlier_weights,
                                estimation, nchunks, progress) {
  
  if (is.null(manifold$library_hrfs)) {
    manifold$library_hrfs <- manifold$B_reconstructor_matrix
  }
  V <- ncol(Y_data)
  chunk_size <- ceiling(V / nchunks)
  
  # Initialize result storage
  Xi_smoothed <- matrix(0, manifold$m_manifold_dim, V)
  H_shapes <- matrix(0, nrow(manifold$B_reconstructor_matrix), V)
  H_raw <- matrix(0, nrow(manifold$B_reconstructor_matrix), V)
  Beta_condition <- matrix(0, design_info$n_conditions, V)
  Beta_condition_initial <- matrix(0, design_info$n_conditions, V)
  Beta_trial <- NULL
  if (estimation %in% c("trial", "both")) {
    Beta_trial <- matrix(0, design_info$n_trials, V)
  }
  
  # Process in chunks
  if (progress) {
    pb <- txtProgressBar(min = 0, max = nchunks, style = 3)
  }
  
  for (chunk in 1:nchunks) {
    # Define chunk indices
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, V)
    chunk_indices <- start_idx:end_idx
    
    # Extract chunk data
    Y_chunk <- Y_data[, chunk_indices, drop = FALSE]
    
    # Extract chunk-specific parameters
    voxel_coords_chunk <- if (!is.null(voxel_coords)) {
      voxel_coords[chunk_indices, , drop = FALSE]
    } else NULL
    
    outlier_weights_chunk <- if (!is.null(outlier_weights)) {
      outlier_weights[, chunk_indices, drop = FALSE]
    } else NULL
    
    # Run standard pipeline on chunk
    chunk_result <- run_mhrf_lss_standard(
      Y_data = Y_chunk,
      design_info = design_info,
      manifold = manifold,
      Z_confounds = Z_confounds,
      voxel_coords = voxel_coords_chunk,
      params = params,
      outlier_weights = outlier_weights_chunk,
      estimation = estimation,
      progress = FALSE  # No nested progress bars
    )
    
    # Store chunk results
    Xi_smoothed[, chunk_indices] <- chunk_result$Xi_smoothed
    H_shapes[, chunk_indices] <- chunk_result$H_shapes
    H_raw[, chunk_indices] <- chunk_result$H_raw
    Beta_condition[, chunk_indices] <- chunk_result$Beta_condition
    Beta_condition_initial[, chunk_indices] <- chunk_result$Beta_condition_initial
    
    if (!is.null(chunk_result$Beta_trial)) {
      Beta_trial[, chunk_indices] <- chunk_result$Beta_trial
    }
    
    if (progress) {
      setTxtProgressBar(pb, chunk)
    }
  }
  
  if (progress) {
    close(pb)
  }
  
  # Aggregate diagnostics
  diagnostics <- list(
    n_voxels_processed = V,
    n_chunks = nchunks,
    chunk_size = chunk_size,
    timestamp = Sys.time()
  )
  
  return(list(
    H_raw = H_raw,
    H_shapes = H_shapes,
    Xi_smoothed = Xi_smoothed,
    Beta_condition_initial = Beta_condition_initial,
    Beta_condition = Beta_condition,
    Beta_trial = Beta_trial,
    diagnostics = diagnostics
  ))
}





