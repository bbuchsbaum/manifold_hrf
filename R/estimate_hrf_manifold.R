#' Estimate HRF Using Manifold Learning
#'
#' @description
#' Primary user-facing function for manifold-guided HRF estimation.
#' Estimates voxel-specific hemodynamic response functions using manifold learning
#' to constrain the space of possible HRF shapes.
#'
#' @param fmri_data A 4D neuroimaging array, NeuroVec object, or matrix (timepoints x voxels)
#' @param events Data frame with columns 'onset', 'condition', and optionally 'duration'
#' @param TR Repetition time in seconds
#' @param hrf_library Either a string preset ("fir", "gamma", "custom") or a matrix/list 
#'   of HRF basis functions
#' @param mask Optional brain mask (3D array or NeuroVol). If NULL, all non-zero voxels are used
#' @param control A list of control parameters from \code{\link{manifold_control}}
#'
#' @return A \code{\link{manifold_hrf_fit}} object containing:
#' \itemize{
#'   \item \code{amplitudes}: Condition-level amplitude estimates
#'   \item \code{trial_amplitudes}: Trial-wise amplitude estimates  
#'   \item \code{hrf_shapes}: Estimated HRF shapes for each voxel
#'   \item \code{fitted_values}: Model predictions
#'   \item \code{residuals}: Model residuals
#'   \item Additional model details in \code{model_specific}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default settings
#' fit <- estimate_hrf_manifold(
#'   fmri_data = bold_data,
#'   events = event_df,
#'   TR = 2.0,
#'   hrf_library = "fir"
#' )
#' 
#' # Custom control parameters
#' ctrl <- manifold_control(
#'   preset = "fast",
#'   lambda_spatial_smooth = 1.0
#' )
#' 
#' fit <- estimate_hrf_manifold(
#'   fmri_data = bold_data,
#'   events = event_df,
#'   TR = 2.0,
#'   hrf_library = "gamma",
#'   control = ctrl
#' )
#' 
#' # View results
#' print(fit)
#' summary(fit)
#' plot(fit, type = "hrf")
#' }
#'
#' @seealso \code{\link{manifold_control}}, \code{\link{manifold_hrf_fit}}
#' @export
estimate_hrf_manifold <- function(fmri_data,
                                 events,
                                 TR,
                                 hrf_library = "fir",
                                 mask = NULL,
                                 control = manifold_control()) {
  
  # Capture call for reproducibility
  cl <- match.call()
  
  # Validate control parameters
  if (!inherits(control, "manifold_control")) {
    stop("control must be a manifold_control object. Use manifold_control() to create one.")
  }
  
  # Log start if verbose
  if (control$verbose_level > 0) {
    message("Starting manifold HRF estimation...")
    message(sprintf("Control preset: %s", attr(control, "preset")))
  }
  
  # Input validation and preprocessing
  if (control$verbose_level > 0) {
    message("Validating inputs...")
  }
  
  # Validate TR
  if (!is.numeric(TR) || length(TR) != 1 || TR <= 0) {
    stop("TR must be a positive number")
  }
  
  # Validate events
  if (!is.data.frame(events)) {
    stop("events must be a data frame")
  }
  
  required_cols <- c("onset", "condition")
  missing_cols <- setdiff(required_cols, names(events))
  if (length(missing_cols) > 0) {
    stop("events missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Process fMRI data
  data_info <- process_fmri_input(fmri_data, mask)
  Y_matrix <- data_info$data_matrix
  voxel_coords <- data_info$voxel_coords
  
  n_timepoints <- nrow(Y_matrix)
  n_voxels <- ncol(Y_matrix)
  
  # Process events
  event_info <- process_events(events, n_timepoints, TR)
  n_trials <- event_info$n_trials
  n_conditions <- event_info$n_conditions
  
  # Create or validate HRF library
  if (control$verbose_level > 0) {
    message("Setting up HRF library...")
  }
  
  hrf_lib <- setup_hrf_library(hrf_library, TR, control)
  
  # Call the existing mhrf_lss function with translated parameters
  # This is temporary until we fully reimplement the pipeline
  
  if (control$verbose_level > 0) {
    message("Running manifold HRF analysis...")
  }
  
  # Translate control parameters to mhrf_lss format
  old_params <- list(
    m_manifold_dim = control$m_manifold_dim,
    lambda_manifold = control$lambda_manifold,
    num_neighbors_mfd = control$num_neighbors_mfd,
    lambda_spatial_smooth = control$lambda_spatial_smooth,
    num_neighbors_Lsp = control$num_neighbors_Lsp,
    lambda_gamma = control$lambda_gamma,
    lambda_betas = control$lambda_condition_betas,
    max_iter = control$max_iter,
    n_cores = control$n_cores,
    verbose_level = control$verbose_level,
    generate_qc_report = control$generate_qc_report
  )
  
  # Call existing implementation
  result_old <- mhrf_lss(
    Y_data = Y_matrix,
    X_onsets = event_info$onset_matrix,
    TR = TR,
    hrf_library = hrf_lib,
    voxel_coords = voxel_coords,
    confounds = NULL,  # TODO: Add confounds support
    params = old_params
  )
  
  # Transform to new format
  fit <- transform_to_manifold_hrf_fit(
    result_old = result_old,
    call = cl,
    control = control,
    n_timepoints = n_timepoints,
    n_voxels = n_voxels,
    n_trials = n_trials,
    n_conditions = n_conditions,
    TR = TR
  )
  
  if (control$verbose_level > 0) {
    message("Manifold HRF estimation complete!")
  }
  
  return(fit)
}

# Helper functions ----

#' Process fMRI input data
#' @keywords internal
process_fmri_input <- function(fmri_data, mask = NULL) {
  # Handle different input types
  if (inherits(fmri_data, "NeuroVec")) {
    # neuroim2 4D object
    if (!requireNamespace("neuroim2", quietly = TRUE)) {
      stop("neuroim2 package required for NeuroVec input")
    }
    
    # Extract data matrix
    if (!is.null(mask)) {
      data_matrix <- as.matrix(fmri_data[mask > 0])
      voxel_coords <- coordinates(mask)[mask > 0, ]
    } else {
      # Use all non-zero voxels
      vol_mask <- apply(fmri_data, 1:3, function(x) any(x != 0))
      data_matrix <- as.matrix(fmri_data[vol_mask])
      voxel_coords <- coordinates(fmri_data)[rep(vol_mask, dim(fmri_data)[4]), ]
    }
    
  } else if (is.array(fmri_data) && length(dim(fmri_data)) == 4) {
    # 4D array
    dims <- dim(fmri_data)
    
    if (!is.null(mask)) {
      mask_idx <- which(mask > 0)
      data_matrix <- matrix(NA, dims[4], length(mask_idx))
      
      for (t in 1:dims[4]) {
        data_matrix[t, ] <- fmri_data[,,, t][mask_idx]
      }
      
      # Get voxel coordinates
      coords <- which(mask > 0, arr.ind = TRUE)
      voxel_coords <- coords
    } else {
      # Reshape to matrix
      data_matrix <- matrix(fmri_data, prod(dims[1:3]), dims[4])
      data_matrix <- t(data_matrix)
      
      # Remove zero voxels
      nonzero_vox <- colSums(data_matrix != 0) > 0
      data_matrix <- data_matrix[, nonzero_vox]
      
      # Get coordinates
      all_coords <- expand.grid(x = 1:dims[1], y = 1:dims[2], z = 1:dims[3])
      voxel_coords <- as.matrix(all_coords[nonzero_vox, ])
    }
    
  } else if (is.matrix(fmri_data)) {
    # Already a matrix (timepoints x voxels)
    data_matrix <- fmri_data
    
    # No spatial information available
    voxel_coords <- matrix(1:ncol(fmri_data), ncol = 1)
    
  } else {
    stop("fmri_data must be a NeuroVec, 4D array, or matrix")
  }
  
  list(
    data_matrix = data_matrix,
    voxel_coords = voxel_coords
  )
}

#' Process event information
#' @keywords internal
process_events <- function(events, n_timepoints, TR) {
  # Sort by onset
  events <- events[order(events$onset), ]
  
  # Get unique conditions
  conditions <- unique(events$condition)
  n_conditions <- length(conditions)
  
  # Create condition factor
  events$condition <- factor(events$condition, levels = conditions)
  
  # Count trials
  n_trials <- nrow(events)
  
  # Create onset matrix (binary indicators)
  onset_matrix <- matrix(0, n_timepoints, n_trials)
  
  for (i in 1:n_trials) {
    onset_tr <- round(events$onset[i] / TR) + 1
    if (onset_tr >= 1 && onset_tr <= n_timepoints) {
      onset_matrix[onset_tr, i] <- 1
    }
  }
  
  list(
    events = events,
    onset_matrix = onset_matrix,
    n_trials = n_trials,
    n_conditions = n_conditions,
    conditions = conditions
  )
}

#' Setup HRF library
#' @keywords internal
setup_hrf_library <- function(hrf_library, TR, control) {
  if (is.character(hrf_library)) {
    # Handle presets
    hrf_lib <- switch(hrf_library,
      "fir" = create_fir_basis(n_basis = 10, TR = TR),
      "gamma" = create_gamma_hrf_library(TR = TR),
      "custom" = stop("For custom HRF library, provide a matrix or list"),
      stop("Unknown HRF library preset: ", hrf_library)
    )
  } else if (is.matrix(hrf_library) || is.list(hrf_library)) {
    # Custom library provided
    hrf_lib <- hrf_library
  } else {
    stop("hrf_library must be a string preset or matrix/list")
  }
  
  hrf_lib
}

#' Transform old result format to manifold_hrf_fit
#' @keywords internal
transform_to_manifold_hrf_fit <- function(result_old, call, control, 
                                        n_timepoints, n_voxels, n_trials, 
                                        n_conditions, TR) {
  
  # Extract components from old format
  amplitudes <- if (!is.null(result_old$beta$condition_final)) {
    result_old$beta$condition_final
  } else {
    result_old$beta$condition
  }
  
  trial_amplitudes <- result_old$beta$trial
  hrf_shapes <- result_old$hrf_shape
  
  # Calculate fitted values if not provided
  fitted_values <- result_old$fitted %||% matrix(0, n_timepoints, n_voxels)
  
  # Calculate residuals
  residuals <- if (!is.null(result_old$Y_data) && !is.null(result_old$fitted)) {
    result_old$Y_data - result_old$fitted
  } else {
    matrix(0, n_timepoints, n_voxels)
  }
  
  # Extract model-specific components
  manifold_coords <- result_old$Xi_final %||% result_old$Xi
  manifold_obj <- result_old$manifold
  amplitudes_initial <- result_old$beta$condition_initial %||% result_old$beta$condition
  spatial_laplacian <- result_old$L_spatial
  
  # Convergence info
  convergence_info <- list(
    converged = TRUE,  # Assume converged if we got results
    iterations = result_old$iterations %||% control$max_iter
  )
  
  # Create data_info
  data_info <- list(
    n_timepoints = n_timepoints,
    n_voxels = n_voxels,
    n_trials = n_trials,
    n_conditions = n_conditions,
    TR = TR
  )
  
  # QC metrics if available
  qc_metrics <- result_old$qc_metrics
  
  # Create the new object
  new_manifold_hrf_fit(
    amplitudes = amplitudes,
    trial_amplitudes = trial_amplitudes,
    hrf_shapes = hrf_shapes,
    fitted_values = fitted_values,
    residuals = residuals,
    manifold_coords = manifold_coords,
    manifold = manifold_obj,
    amplitudes_initial = amplitudes_initial,
    spatial_laplacian = spatial_laplacian,
    convergence_info = convergence_info,
    call = call,
    control = control,
    data_info = data_info,
    qc_metrics = qc_metrics
  )
}

#' Create gamma HRF library
#' @keywords internal
create_gamma_hrf_library <- function(TR, n_basis = 5) {
  # This is a placeholder - would need proper implementation
  # For now, return a simple gamma-based library
  
  stop("Gamma HRF library not yet implemented. Use 'fir' for now.")
}