#' Simple Parameter Helper for mhrf_lss
#'
#' Creates a parameter list for `mhrf_lss` using preset defaults
#' from \code{get_preset_params()} and allowing user overrides.
#'
#' @param preset Character preset passed to \code{get_preset_params}.
#' @param ... Additional parameter overrides.
#' @return A named list of parameters.
#' @export
mhrf_lss_parameters <- function(preset = "balanced", ...) {
  params <- get_preset_params(preset)
  user <- list(...)
  for (nm in names(user)) {
    params[[nm]] <- user[[nm]]
  }
  params
}

#' Run the Core M-HRF-LSS Pipeline
#'
#' This wrapper provides a streamlined interface for running the
#' M-HRF-LSS algorithm when design matrices are already available.
#'
#' @param Y_bold Numeric matrix of BOLD time series (time \eqn{\times} voxels).
#' @param X_condition_list List of condition-level design matrices.
#' @param X_trial_list List of trial-level design matrices.
#' @param Z_confounds Optional matrix of confound regressors.
#' @param voxel_coordinates Optional matrix of voxel coordinates for
#'   spatial smoothing.
#' @param TR Repetition time in seconds.
#' @param parameters List of algorithm parameters as created by
#'   \code{mhrf_lss_parameters()}.
#' @return List with M-HRF-LSS results.
#' @export
mhrf_lss <- function(Y_bold,
                     X_condition_list,
                     X_trial_list,
                     Z_confounds = NULL,
                     voxel_coordinates = NULL,
                     TR = 2,
                     parameters = mhrf_lss_parameters()) {

  stopifnot(is.matrix(Y_bold))
  stopifnot(is.list(X_condition_list))
  stopifnot(is.list(X_trial_list))

  design_info <- list(
    X_condition_list = X_condition_list,
    X_trial_list = X_trial_list,
    n_conditions = length(X_condition_list),
    n_trials = length(X_trial_list)
  )

  manifold <- create_hrf_manifold(
    hrf_library = "gamma_grid",
    params = parameters,
    TR = TR,
    verbose = FALSE
  )

  result <- run_mhrf_lss_standard(
    Y_data = Y_bold,
    design_info = design_info,
    manifold = manifold,
    Z_confounds = Z_confounds,
    voxel_coords = voxel_coordinates,
    params = parameters,
    outlier_weights = NULL,
    estimation = if (length(X_trial_list) > 0) "both" else "condition",
    progress = FALSE
  )

  return(result)
}

