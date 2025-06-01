#' Plot QC summary for mhrf_result
#'
#' Generates the default QC summary plot used by `plot()` for
#' `mhrf_result` objects, but as a standalone function.
#'
#' @param x mhrf_result object
#' @return Invisibly returns `x`.
#' @export
plot_qc_summary <- function(x) {
  stopifnot(inherits(x, "mhrf_result"))
  .plot_quality_metrics(x)
  invisible(x)
}
