#' Simulate a Simple fMRI Dataset
#'
#' Generates synthetic BOLD data and accompanying design information
#' using the internal simulation helpers. The output is wrapped in a
#' `fmrireg::matrix_dataset` object so it can be used directly with
#' functions that expect `fmrireg` datasets.
#'
#' @param n_voxels Number of voxels to simulate.
#' @param n_timepoints Number of time points.
#' @param n_trials Number of trials per condition.
#' @param n_conditions Number of experimental conditions.
#' @param TR Repetition time in seconds.
#' @param hrf_variability HRF variability level: "none", "moderate", or "high".
#' @param noise_level DVARS noise level as percentage.
#' @param seed Random seed for reproducibility.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{dataset}: `fmrireg` matrix_dataset with the simulated BOLD data
#'     \item \code{bold_data}: Matrix of noisy BOLD signals (time x voxel)
#'     \item \code{confounds}: Matrix of confound regressors
#'     \item \code{ground_truth}: List with HRFs, amplitudes and design info
#'   }
#'
#' @examples
#' \dontrun{
#' sim <- simulate_mhrf_dataset(n_voxels = 50, n_timepoints = 150)
#' }
#'
#' @export
simulate_mhrf_dataset <- function(n_voxels = 50,
                                  n_timepoints = 150,
                                  n_trials = 10,
                                  n_conditions = 2,
                                  TR = 1.0,
                                  hrf_variability = c("none", "moderate", "high"),
                                  noise_level = 5,
                                  seed = 1) {

  set.seed(seed)
  hrf_variability <- match.arg(hrf_variability)

  # Ground truth HRFs
  hrfs <- generate_ground_truth_hrfs(
    n_voxels = n_voxels,
    hrf_variability = hrf_variability,
    TR = TR,
    manifold_params = list(TR_precision = 0.1)
  )

  # Experimental design
  design <- generate_experimental_design(
    n_timepoints = n_timepoints,
    n_trials = n_trials,
    n_conditions = n_conditions,
    TR = TR,
    hrf_length = length(hrfs$time_points)
  )

  # Amplitudes
  amps <- generate_ground_truth_amplitudes(
    n_voxels = n_voxels,
    n_conditions = n_conditions,
    n_trials = design$total_trials,
    activation_patterns = c("sustained", "transient", "mixed")
  )

  # BOLD data with noise
  bold <- generate_bold_data(
    ground_truth_hrfs = hrfs,
    ground_truth_amplitudes = amps,
    design_info = design,
    noise_level = noise_level,
    TR = TR
  )

  # Event table for fmrireg dataset
  event_table <- data.frame(
    onset = design$onsets,
    duration = 0,
    amplitude = 1,
    condition = factor(design$conditions),
    trial = seq_along(design$onsets),
    run = 1
  )

  dset <- fmrireg::matrix_dataset(
    datamat = bold$Y_data,
    TR = TR,
    run_length = n_timepoints,
    event_table = event_table
  )

  list(
    dataset = dset,
    bold_data = bold$Y_data,
    confounds = bold$Z_confounds,
    ground_truth = list(
      hrfs = hrfs,
      amplitudes = amps,
      design = design
    )
  )
}
