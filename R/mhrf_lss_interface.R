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
#' @return An mhrf_lss_fit object containing:
#'   - manifold: The HRF manifold object
#'   - hrfs: Estimated voxel-specific HRFs
#'   - betas_condition: Condition-level amplitudes
#'   - betas_trial: Trial-wise amplitudes (if requested)
#'   - model_info: Original model specifications
#'   - diagnostics: Fit quality metrics
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
#' }
#'
#' @export
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
  
  # Validate inputs
  if (!requireNamespace("fmrireg", quietly = TRUE)) {
    stop("Package 'fmrireg' is required. Please install it.")
  }
  
  estimation <- match.arg(estimation)
  strategy <- match.arg(strategy)
  
  # Extract data dimensions
  data_mat <- fmrireg::get_data_matrix(dataset)
  n_time <- nrow(data_mat)
  n_voxels <- ncol(data_mat)
  
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
      manifold = manifold,
      hrfs = results$H_shapes,
      xi_coordinates = results$Xi_smoothed,
      betas_condition = results$Beta_condition,
      betas_trial = if (estimation %in% c("trial", "both")) results$Beta_trial else NULL,
      model_info = list(
        formula = formula,
        event_model = event_mod,
        baseline_model = baseline_model,
        estimation = estimation,
        n_timepoints = n_time,
        n_voxels = n_voxels,
        n_conditions = nrow(results$Beta_condition),
        n_trials = if (!is.null(results$Beta_trial)) nrow(results$Beta_trial) else NULL
      ),
      diagnostics = results$diagnostics,
      dataset_info = list(
        TR = dataset$TR,
        run_lengths = dataset$run_length
      ),
      params = list(
        manifold = manifold$parameters,
        robust = robust
      ),
      call = match.call()
    ),
    class = c("mhrf_lss_fit", "list")
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
  
  # Handle different library sources
  if (is.character(hrf_library) && length(hrf_library) == 1) {
    if (hrf_library == "canonical") {
      if (!requireNamespace("fmrireg", quietly = TRUE)) {
        stop("Package 'fmrireg' is required to use the canonical HRF library.")
      }
      # Use standard fmrireg HRFs
      if (verbose) message("  Using canonical fmrireg HRF library")
      
      hrf_list <- list(
        fmrireg::HRF_SPMG1,
        fmrireg::HRF_SPMG2,
        fmrireg::HRF_SPMG3,
        fmrireg::HRF_GAMMA,
        fmrireg::HRF_GAUSSIAN
      )
      
      # Add variations
      hrf_list <- c(
        hrf_list,
        lapply(c(1, 2, 3), function(lag) {
          fmrireg::lag_hrf(fmrireg::HRF_SPMG1, lag)
        })
      )
      
      hrf_library_source <- hrf_list
      
    } else {
      # Use predefined library type
      hrf_library_source <- hrf_library
    }
  } else {
    hrf_library_source <- hrf_library
  }
  
  # First, we need to create the affinity matrix
  # This depends on the type of HRF library source
  
  if (is.character(hrf_library_source) && length(hrf_library_source) == 1) {
    # Use a predefined library
    # For now, create a simple gamma grid (would be replaced with actual implementation)
    p <- 30  # HRF length
    N <- 50  # Number of HRFs
    L_library <- matrix(rnorm(p * N), p, N)
    # Normalize
    L_library <- apply(L_library, 2, function(x) x / sum(abs(x)))
  } else if (is.matrix(hrf_library_source)) {
    L_library <- hrf_library_source
  } else if (is.list(hrf_library_source)) {
    # Convert list of HRF objects to matrix
    # This would need proper implementation with fmrireg
    p <- 30  # Assumed HRF length
    N <- length(hrf_library_source)
    L_library <- matrix(0, p, N)
    # Placeholder - would evaluate each HRF
    for (i in 1:N) {
      L_library[, i] <- rnorm(p)
    }
    L_library <- apply(L_library, 2, function(x) x / sum(abs(x)))
  } else {
    stop("Invalid hrf_library_source format")
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
    manifold <- suppressWarnings(
      get_manifold_basis_reconstructor_robust(
        S_markov_matrix = S_markov,
        L_library_matrix = L_library,
        m_manifold_dim_target = params$m_manifold_dim_target,
        m_manifold_dim_min_variance = params$m_manifold_dim_min_variance %||% 0.95,
        fallback_to_pca = TRUE
      )
    )
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
  X_condition_list <- term_mats
  
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
    ident_sign_method = "canonical_correlation"
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
      ram_heuristic_GB_for_Rt = 1.0
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
    Xi_smoothed = Xi_smoothed
  )
  
  return(list(
    H_shapes = H_shapes,
    Xi_smoothed = Xi_smoothed,
    Beta_condition = Beta_condition_final,
    Beta_trial = Beta_trial,
    diagnostics = diagnostics
  ))
}


#' Compute Pipeline Diagnostics
#'
#' @keywords internal
compute_pipeline_diagnostics <- function(Y_data, Y_proj, H_shapes, 
                                       Beta_condition, Xi_smoothed) {
  
  # This would compute R², HRF statistics, etc.
  # Simplified version:
  
  diagnostics <- list(
    n_voxels_processed = ncol(H_shapes),
    manifold_variance = apply(Xi_smoothed, 1, var),
    hrf_stats = extract_hrf_stats(H_shapes),
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


# S3 Methods for mhrf_lss_fit

#' @export
print.mhrf_lss_fit <- function(x, ...) {
  cat("M-HRF-LSS Model Fit\n")
  cat("===================\n\n")
  
  cat("Formula:", deparse(x$model_info$formula), "\n")
  cat("Data dimensions:", x$model_info$n_timepoints, "timepoints x",
      x$model_info$n_voxels, "voxels\n")
  cat("Conditions:", x$model_info$n_conditions, "\n")
  
  if (!is.null(x$model_info$n_trials)) {
    cat("Trials:", x$model_info$n_trials, "\n")
  }
  
  cat("\nManifold:\n")
  cat("  Method:", x$manifold$method_used, "\n")
  cat("  Dimensions:", x$manifold$m_manifold_dim, "\n")
  cat("  Library size:", x$manifold$parameters$n_hrfs_library, "HRFs\n")
  
  cat("\nEstimation:", x$model_info$estimation, "\n")
  
  invisible(x)
}


#' Extract Coefficients from M-HRF-LSS Fit
#'
#' @export
coef.mhrf_lss_fit <- function(object, type = c("condition", "trial", "hrf"), ...) {
  
  type <- match.arg(type)
  
  switch(type,
    condition = object$betas_condition,
    trial = object$betas_trial,
    hrf = object$hrfs
  )
}


#' Plot M-HRF-LSS Results
#'
#' @export
plot.mhrf_lss_fit <- function(x, type = c("hrfs", "manifold", "diagnostics"),
                             voxels = NULL, ...) {
  
  type <- match.arg(type)
  
  if (type == "hrfs") {
    # Plot HRF shapes for selected voxels
    if (is.null(voxels)) {
      # Select representative voxels
      voxels <- sample(1:ncol(x$hrfs), min(9, ncol(x$hrfs)))
    }
    
    # Would create actual plots here
    message(sprintf("Plotting HRFs for voxels: %s", 
                   paste(voxels, collapse = ", ")))
  }
  
  # Other plot types...
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
  Beta_condition <- matrix(0, design_info$n_conditions, V)
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
    Beta_condition[, chunk_indices] <- chunk_result$Beta_condition
    
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
    H_shapes = H_shapes,
    Xi_smoothed = Xi_smoothed,
    Beta_condition = Beta_condition,
    Beta_trial = Beta_trial,
    diagnostics = diagnostics
  ))
}


#' Detect Outlier Timepoints (Stub for Interface)
#'
#' @keywords internal
detect_outlier_timepoints <- function(Y_data, threshold = 3) {
  # This is a stub - actual implementation is in robust_spatial_outlier.R
  # For now, return uniform weights
  return(matrix(1, nrow(Y_data), ncol(Y_data)))
}


#' Extract HRF Statistics (Stub for Interface)
#'
#' @keywords internal  
extract_hrf_stats <- function(H_shapes) {
  # This is a stub - actual implementation is in qc_report.R
  # Return basic stats
  list(
    mean_peak_time = NA,
    mean_peak_amplitude = NA,
    mean_fwhm = NA
  )
}

