#' Plot Trial-wise Amplitudes for a Voxel
#'
#' Visualize estimated trial-level amplitudes for a single voxel from an
#' `mhrf_result` object. Optionally overlay the ground-truth amplitudes
#' for comparison.
#'
#' @param result An object of class `mhrf_result` containing trial estimates.
#' @param voxel Integer index of the voxel to plot.
#' @param true_values Optional numeric vector of true trial amplitudes to
#'   overlay on the plot.
#' @param ... Additional graphical parameters passed to [graphics::plot].
#'
#' @return Invisibly returns the vector of trial amplitudes for the selected
#'   voxel.
#' @export
plot_trial_betas <- function(result, voxel = 1, true_values = NULL, ...) {
  if (!inherits(result, "mhrf_result")) {
    stop("'result' must be an mhrf_result object")
  }

  betas <- coef(result, type = "trial_amplitudes")
  if (is.null(betas)) {
    stop("Trial amplitudes not available in result")
  }

  if (voxel < 1 || voxel > ncol(betas)) {
    stop("voxel index out of range")
  }

  est <- betas[, voxel]
  n_trials <- length(est)

  plot(seq_len(n_trials), est, type = "b", pch = 16, col = "red",
       xlab = "Trial", ylab = "Amplitude",
       main = sprintf("Trial-wise Amplitudes (voxel %d)", voxel), ...)
  abline(h = 0, lty = 2, col = "gray")

  if (!is.null(true_values)) {
    if (length(true_values) != n_trials) {
      stop("true_values must have length equal to number of trials")
    }
    lines(seq_len(n_trials), true_values, type = "b", pch = 1,
          col = "blue", lty = 2)
    legend("topright", legend = c("Estimated", "True"),
           col = c("red", "blue"), pch = c(16, 1), lty = c(1, 2))
  } else {
    legend("topright", legend = "Estimated", col = "red", pch = 16, lty = 1)
  }

  invisible(est)
}
