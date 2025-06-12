#' S3 Methods for manifold_hrf_fit Objects
#'
#' @description
#' A collection of S3 methods for working with manifold_hrf_fit objects,
#' including print, summary, coef, plot, and prediction methods.
#'
#' @name manifold_hrf_fit_methods
NULL

#' Print method for manifold_hrf_fit objects
#'
#' @param x A manifold_hrf_fit object
#' @param digits Number of digits to display
#' @param ... Additional arguments (ignored)
#' @export
print.manifold_hrf_fit <- function(x, digits = 3, ...) {
  cat("Manifold HRF Fit\n")
  cat("================\n")
  
  # Call
  cat("\nCall:\n")
  print(x$call)
  
  # Data dimensions
  cat("\nData:\n")
  cat(sprintf("  Timepoints: %d\n", x$data_info$n_timepoints))
  cat(sprintf("  Voxels: %d\n", x$data_info$n_voxels))
  cat(sprintf("  Trials: %d\n", x$data_info$n_trials))
  cat(sprintf("  Conditions: %d\n", x$data_info$n_conditions))
  cat(sprintf("  TR: %.2f seconds\n", x$data_info$TR))
  
  # Model info
  cat("\nModel:\n")
  cat(sprintf("  Manifold dimension: %d\n", ncol(x$model_specific$manifold_coords)))
  cat(sprintf("  HRF duration: %.1f seconds (%d samples)\n", 
              nrow(x$hrf_shapes) * x$data_info$TR, nrow(x$hrf_shapes)))
  
  # Convergence
  if (!is.null(x$model_specific$convergence_info)) {
    conv <- x$model_specific$convergence_info
    cat(sprintf("  Converged: %s", ifelse(conv$converged, "Yes", "No")))
    if (!is.null(conv$iterations)) {
      cat(sprintf(" (iterations: %d)", conv$iterations))
    }
    cat("\n")
  }
  
  # Condition amplitudes summary
  cat("\nCondition Amplitudes:\n")
  if (is.matrix(x$amplitudes)) {
    # Multi-voxel case
    amp_summary <- apply(x$amplitudes, 1, summary)
    rownames(amp_summary) <- c("Min", "Q1", "Median", "Mean", "Q3", "Max")
    print(round(amp_summary, digits))
  } else {
    # Single voxel case
    print(round(x$amplitudes, digits))
  }
  
  # Model fit
  if (!is.null(x$fitted_values) && !is.null(x$residuals)) {
    r_squared <- 1 - sum(x$residuals^2) / sum((x$fitted_values + x$residuals - mean(x$fitted_values + x$residuals))^2)
    cat(sprintf("\nModel fit: R² = %.3f\n", r_squared))
  }
  
  invisible(x)
}

#' Summary method for manifold_hrf_fit objects
#'
#' @param object A manifold_hrf_fit object
#' @param ... Additional arguments (ignored)
#' @return A summary.manifold_hrf_fit object
#' @export
summary.manifold_hrf_fit <- function(object, ...) {
  # Calculate summary statistics
  r_squared <- if (!is.null(object$fitted_values) && !is.null(object$residuals)) {
    1 - sum(object$residuals^2) / sum((object$fitted_values + object$residuals - mean(object$fitted_values + object$residuals))^2)
  } else {
    NA
  }
  
  # Amplitude summaries
  if (is.matrix(object$amplitudes)) {
    amplitude_summary <- t(apply(object$amplitudes, 1, summary))
  } else {
    amplitude_summary <- summary(object$amplitudes)
  }
  
  # Trial amplitude summary
  if (is.matrix(object$trial_amplitudes)) {
    trial_amplitude_summary <- summary(as.vector(object$trial_amplitudes))
  } else {
    trial_amplitude_summary <- summary(object$trial_amplitudes)
  }
  
  # HRF characteristics
  hrf_peaks <- apply(object$hrf_shapes, 2, which.max) * object$data_info$TR
  hrf_peak_summary <- summary(hrf_peaks)
  
  structure(
    list(
      call = object$call,
      data_info = object$data_info,
      control_preset = attr(object$control, "preset"),
      manifold_dim = ncol(object$model_specific$manifold_coords),
      amplitude_summary = amplitude_summary,
      trial_amplitude_summary = trial_amplitude_summary,
      hrf_peak_summary = hrf_peak_summary,
      convergence = object$model_specific$convergence_info,
      r_squared = r_squared,
      qc_available = !is.null(object$qc_metrics)
    ),
    class = "summary.manifold_hrf_fit"
  )
}

#' Print method for summary.manifold_hrf_fit
#'
#' @param x A summary.manifold_hrf_fit object
#' @param digits Number of digits to display
#' @param ... Additional arguments (ignored)
#' @export
print.summary.manifold_hrf_fit <- function(x, digits = 3, ...) {
  cat("Summary of Manifold HRF Fit\n")
  cat("===========================\n")
  
  cat("\nCall:\n")
  print(x$call)
  
  cat("\nData Dimensions:\n")
  cat(sprintf("  Timepoints: %d\n", x$data_info$n_timepoints))
  cat(sprintf("  Voxels: %d\n", x$data_info$n_voxels))
  cat(sprintf("  Trials: %d\n", x$data_info$n_trials))
  cat(sprintf("  Conditions: %d\n", x$data_info$n_conditions))
  
  cat("\nModel Configuration:\n")
  cat(sprintf("  Control preset: %s\n", x$control_preset))
  cat(sprintf("  Manifold dimension: %d\n", x$manifold_dim))
  
  cat("\nCondition Amplitude Summary:\n")
  print(round(x$amplitude_summary, digits))
  
  cat("\nTrial Amplitude Summary:\n")
  print(round(x$trial_amplitude_summary, digits))
  
  cat("\nHRF Peak Time Summary (seconds):\n")
  print(round(x$hrf_peak_summary, digits))
  
  if (!is.null(x$convergence)) {
    cat("\nConvergence:\n")
    cat(sprintf("  Status: %s\n", ifelse(x$convergence$converged, "Converged", "Not converged")))
    if (!is.null(x$convergence$iterations)) {
      cat(sprintf("  Iterations: %d\n", x$convergence$iterations))
    }
  }
  
  if (!is.na(x$r_squared)) {
    cat(sprintf("\nModel Fit: R² = %.3f\n", x$r_squared))
  }
  
  if (x$qc_available) {
    cat("\nQuality control metrics available - use qc_report() to view\n")
  }
  
  invisible(x)
}

#' Extract coefficients from manifold_hrf_fit
#'
#' @param object A manifold_hrf_fit object
#' @param type Type of coefficients to extract
#' @param ... Additional arguments (ignored)
#' @return The requested coefficients
#' @export
coef.manifold_hrf_fit <- function(object, 
                                 type = c("amplitudes", "trial_amplitudes", "manifold_coords"),
                                 ...) {
  type <- match.arg(type)
  
  switch(type,
    amplitudes = object$amplitudes,
    trial_amplitudes = object$trial_amplitudes,
    manifold_coords = object$model_specific$manifold_coords
  )
}

#' Plot manifold_hrf_fit objects
#'
#' @param x A manifold_hrf_fit object
#' @param type Type of plot to create
#' @param voxel For HRF plots, which voxel(s) to plot (default: 1:5)
#' @param conditions For amplitude plots, which conditions to show
#' @param ... Additional plotting arguments
#' @export
plot.manifold_hrf_fit <- function(x, 
                                 type = c("hrf", "manifold", "amplitudes", "fit"),
                                 voxel = NULL,
                                 conditions = NULL,
                                 ...) {
  type <- match.arg(type)
  
  switch(type,
    hrf = plot_hrf_shapes(x, voxel = voxel, ...),
    manifold = plot_manifold_coords(x, ...),
    amplitudes = plot_amplitude_distribution(x, conditions = conditions, ...),
    fit = plot_model_fit(x, voxel = voxel, ...)
  )
}

#' Fitted values from manifold_hrf_fit
#'
#' @param object A manifold_hrf_fit object
#' @param ... Additional arguments (ignored)
#' @return The fitted values
#' @export
fitted.manifold_hrf_fit <- function(object, ...) {
  object$fitted_values
}

#' Residuals from manifold_hrf_fit
#'
#' @param object A manifold_hrf_fit object
#' @param type Type of residuals (default: "response")
#' @param ... Additional arguments (ignored)
#' @return The residuals
#' @export
residuals.manifold_hrf_fit <- function(object, type = "response", ...) {
  if (type != "response") {
    warning("Only 'response' residuals are currently implemented")
  }
  object$residuals
}

#' Predict method for manifold_hrf_fit
#'
#' @param object A manifold_hrf_fit object
#' @param newdata New event data frame (optional)
#' @param type Type of prediction ("response" or "hrf")
#' @param ... Additional arguments
#' @return Predicted values
#' @export
predict.manifold_hrf_fit <- function(object, newdata = NULL, type = "response", ...) {
  if (is.null(newdata)) {
    # Return fitted values for original data
    return(fitted(object))
  }
  
  # Validate newdata
  if (!is.data.frame(newdata)) {
    stop("newdata must be a data frame")
  }
  
  required_cols <- c("onset", "condition")
  missing_cols <- setdiff(required_cols, names(newdata))
  if (length(missing_cols) > 0) {
    stop("newdata missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Generate predictions for new events
  # This is a placeholder - actual implementation would need to:
  # 1. Create design matrices from new events
  # 2. Convolve with estimated HRFs
  # 3. Apply estimated amplitudes
  
  stop("Prediction for new data not yet implemented")
}

# Internal plotting functions ----

#' Plot HRF shapes
#' @keywords internal
plot_hrf_shapes <- function(x, voxel = NULL, ...) {
  if (is.null(voxel)) {
    # Default to first 5 voxels or all if fewer
    voxel <- seq_len(min(5, x$data_info$n_voxels))
  }
  
  # Extract HRF data
  hrf_data <- x$hrf_shapes[, voxel, drop = FALSE]
  time_points <- (seq_len(nrow(hrf_data)) - 1) * x$data_info$TR
  
  # Create plot
  matplot(time_points, hrf_data, type = "l", lty = 1,
          xlab = "Time (seconds)", ylab = "HRF Amplitude",
          main = "Estimated HRF Shapes", ...)
  abline(h = 0, col = "gray", lty = 2)
  
  if (length(voxel) <= 10) {
    legend("topright", legend = paste("Voxel", voxel), 
           col = seq_along(voxel), lty = 1, bty = "n")
  }
}

#' Plot manifold coordinates
#' @keywords internal
plot_manifold_coords <- function(x, dims = c(1, 2), ...) {
  coords <- x$model_specific$manifold_coords
  
  if (max(dims) > ncol(coords)) {
    stop("Requested dimensions exceed manifold dimension")
  }
  
  plot(coords[dims[1], ], coords[dims[2], ],
       xlab = paste("Manifold Dimension", dims[1]),
       ylab = paste("Manifold Dimension", dims[2]),
       main = "Manifold Coordinates",
       pch = 19, ...)
}

#' Plot amplitude distributions
#' @keywords internal
plot_amplitude_distribution <- function(x, conditions = NULL, ...) {
  if (is.null(conditions)) {
    conditions <- seq_len(x$data_info$n_conditions)
  }
  
  if (is.matrix(x$amplitudes)) {
    # Multi-voxel: show distribution across voxels for each condition
    amp_data <- x$amplitudes[conditions, , drop = FALSE]
    boxplot(t(amp_data), names = paste("Cond", conditions),
            xlab = "Condition", ylab = "Amplitude",
            main = "Amplitude Distribution Across Voxels", ...)
  } else {
    # Single voxel: show bar plot
    barplot(x$amplitudes[conditions], names.arg = paste("Cond", conditions),
            xlab = "Condition", ylab = "Amplitude",
            main = "Condition Amplitudes", ...)
  }
}

#' Plot model fit
#' @keywords internal
plot_model_fit <- function(x, voxel = 1, time_range = NULL, ...) {
  if (is.null(x$fitted_values) || is.null(x$residuals)) {
    stop("Model fit components not available")
  }
  
  # Extract data for specified voxel
  observed <- x$fitted_values[, voxel] + x$residuals[, voxel]
  fitted <- x$fitted_values[, voxel]
  
  # Time range
  if (is.null(time_range)) {
    time_range <- c(1, min(200, length(observed)))
  }
  idx <- time_range[1]:time_range[2]
  
  # Create plot
  plot(idx, observed[idx], type = "l", col = "black",
       xlab = "Time (TRs)", ylab = "BOLD Signal",
       main = paste("Model Fit - Voxel", voxel), ...)
  lines(idx, fitted[idx], col = "red", lwd = 2)
  legend("topright", legend = c("Observed", "Fitted"), 
         col = c("black", "red"), lty = 1, lwd = c(1, 2), bty = "n")
}