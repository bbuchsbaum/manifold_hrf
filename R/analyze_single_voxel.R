#' Analyze a Single Voxel with One Condition
#'
#' Convenience wrapper around `mhrf_analyze` for the simplest possible
#' case of one voxel and one condition. It serves as a small sanity
#' check that the M-HRF-LSS pipeline works end-to-end on minimal data.
#'
#' @param y Numeric vector representing the BOLD time series.
#' @param onsets Numeric vector of event onsets (in seconds).
#' @param TR Repetition time in seconds. Defaults to 2.
#' @param preset Parameter preset passed to `mhrf_analyze`.
#' @param hrf_library HRF library specification for `mhrf_analyze`.
#'
#' @return An object of class `mhrf_result` for the single voxel analysis.
#' @examples
#' \dontrun{
#' y <- rnorm(120)
#' res <- analyze_single_voxel(y, c(20, 60, 100), TR = 1)
#' }
#' @export
analyze_single_voxel <- function(y,
                                 onsets,
                                 TR = 2,
                                 preset = "balanced",
                                 hrf_library = "auto") {
  if (!is.numeric(y)) {
    stop("y must be a numeric vector")
  }
  if (!is.numeric(onsets)) {
    stop("onsets must be a numeric vector")
  }

  events <- data.frame(
    onset = onsets,
    duration = 0,
    condition = "cond1"
  )

  result <- mhrf_analyze(
    Y_data = matrix(y, ncol = 1),
    events = events,
    TR = TR,
    preset = preset,
    hrf_library = hrf_library,
    voxel_mask = TRUE,
    n_jobs = 1,
    verbose = FALSE
  )

  return(result)
}
