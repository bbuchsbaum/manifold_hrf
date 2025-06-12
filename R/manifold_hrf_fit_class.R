#' Manifold HRF Fit Object
#'
#' @description
#' S3 class representing the results of manifold HRF estimation.
#' This object contains the estimated HRF shapes, condition and trial amplitudes,
#' model fit statistics, and technical details about the manifold construction.
#'
#' @details
#' A \code{manifold_hrf_fit} object contains the following components:
#' \describe{
#'   \item{amplitudes}{Condition-level amplitude estimates}
#'   \item{trial_amplitudes}{Trial-wise amplitude estimates}
#'   \item{hrf_shapes}{Matrix of estimated HRF shapes (timepoints x voxels)}
#'   \item{fitted_values}{Predicted BOLD signal}
#'   \item{residuals}{Residuals (observed - fitted)}
#'   \item{model_specific}{List containing technical details:
#'     \itemize{
#'       \item manifold_coords: Manifold coordinates (Xi)
#'       \item manifold: The manifold object
#'       \item amplitudes_initial: Initial amplitude estimates
#'       \item spatial_laplacian: Spatial graph Laplacian
#'       \item convergence_info: Convergence information
#'     }
#'   }
#'   \item{call}{The matched call}
#'   \item{control}{Control parameters used}
#'   \item{data_info}{Data dimensions and metadata}
#'   \item{qc_metrics}{Optional quality control metrics}
#' }
#'
#' @seealso \code{\link{estimate_hrf_manifold}}, \code{\link{manifold_control}}
#' @name manifold_hrf_fit
NULL

#' Create a manifold_hrf_fit object
#'
#' @description
#' Internal function to create a properly structured manifold_hrf_fit object.
#'
#' @param amplitudes Condition-level amplitudes
#' @param trial_amplitudes Trial-wise amplitudes
#' @param hrf_shapes Estimated HRF shapes
#' @param fitted_values Fitted BOLD values
#' @param residuals Model residuals
#' @param manifold_coords Manifold coordinates
#' @param manifold Manifold object
#' @param amplitudes_initial Initial amplitude estimates
#' @param spatial_laplacian Spatial graph Laplacian
#' @param convergence_info Convergence information
#' @param call Matched call
#' @param control Control parameters
#' @param data_info Data information list
#' @param qc_metrics Optional QC metrics
#'
#' @return A manifold_hrf_fit object
#' @keywords internal
new_manifold_hrf_fit <- function(amplitudes,
                                trial_amplitudes,
                                hrf_shapes,
                                fitted_values,
                                residuals,
                                manifold_coords,
                                manifold,
                                amplitudes_initial,
                                spatial_laplacian,
                                convergence_info,
                                call,
                                control,
                                data_info,
                                qc_metrics = NULL) {
  
  # Validate essential components
  if (!is.numeric(amplitudes)) {
    stop("amplitudes must be numeric")
  }
  
  if (!is.matrix(trial_amplitudes) && !is.numeric(trial_amplitudes)) {
    stop("trial_amplitudes must be numeric matrix or vector")
  }
  
  if (!is.matrix(hrf_shapes)) {
    stop("hrf_shapes must be a matrix")
  }
  
  # Create the object
  structure(
    list(
      # Primary results
      amplitudes = amplitudes,
      trial_amplitudes = trial_amplitudes,
      hrf_shapes = hrf_shapes,
      
      # Model fit
      fitted_values = fitted_values,
      residuals = residuals,
      
      # Model details
      model_specific = list(
        manifold_coords = manifold_coords,
        manifold = manifold,
        amplitudes_initial = amplitudes_initial,
        spatial_laplacian = spatial_laplacian,
        convergence_info = convergence_info
      ),
      
      # Metadata
      call = call,
      control = control,
      data_info = data_info,
      
      # Optional
      qc_metrics = qc_metrics
    ),
    class = "manifold_hrf_fit"
  )
}

#' Check if object is a manifold_hrf_fit
#'
#' @param x Object to test
#' @return Logical indicating if x is a manifold_hrf_fit object
#' @export
is.manifold_hrf_fit <- function(x) {
  inherits(x, "manifold_hrf_fit")
}

#' Validate manifold_hrf_fit object
#'
#' @description
#' Internal function to validate the structure of a manifold_hrf_fit object.
#'
#' @param x A manifold_hrf_fit object
#' @return TRUE if valid, otherwise throws an error
#' @keywords internal
validate_manifold_hrf_fit <- function(x) {
  if (!is.manifold_hrf_fit(x)) {
    stop("Object is not a manifold_hrf_fit")
  }
  
  required_components <- c("amplitudes", "trial_amplitudes", "hrf_shapes",
                          "fitted_values", "residuals", "model_specific",
                          "call", "control", "data_info")
  
  missing <- setdiff(required_components, names(x))
  if (length(missing) > 0) {
    stop("Missing required components: ", paste(missing, collapse = ", "))
  }
  
  # Validate model_specific
  required_model_specific <- c("manifold_coords", "manifold", "convergence_info")
  missing_ms <- setdiff(required_model_specific, names(x$model_specific))
  if (length(missing_ms) > 0) {
    stop("Missing required model_specific components: ", 
         paste(missing_ms, collapse = ", "))
  }
  
  # Validate data_info
  required_data_info <- c("n_timepoints", "n_voxels", "n_trials", "n_conditions", "TR")
  missing_di <- setdiff(required_data_info, names(x$data_info))
  if (length(missing_di) > 0) {
    stop("Missing required data_info components: ", 
         paste(missing_di, collapse = ", "))
  }
  
  # Check dimensions consistency
  n_voxels <- x$data_info$n_voxels
  n_conditions <- x$data_info$n_conditions
  n_trials <- x$data_info$n_trials
  
  if (length(x$amplitudes) != n_conditions * n_voxels &&
      !is.matrix(x$amplitudes) || 
      (is.matrix(x$amplitudes) && nrow(x$amplitudes) != n_conditions)) {
    stop("amplitudes dimensions inconsistent with data_info")
  }
  
  if (ncol(x$hrf_shapes) != n_voxels) {
    stop("hrf_shapes columns must equal n_voxels")
  }
  
  TRUE
}