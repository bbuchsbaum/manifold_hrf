#' S3 Methods for mhrf_result Objects
#'
#' Methods for working with results from the M-HRF-LSS analysis

#' Print method for mhrf_result
#'
#' @param x An mhrf_result object
#' @param ... Additional arguments (ignored)
#' @export
print.mhrf_result <- function(x, ...) {
  cat("\nM-HRF-LSS Result\n")
  cat("================\n\n")
  
  # Data summary
  cat("Data:\n")
  cat(sprintf("  • %d timepoints × %d voxels analyzed\n", 
              x$metadata$n_timepoints, x$metadata$n_voxels_analyzed))
  
  if (x$metadata$n_voxels_input != x$metadata$n_voxels_analyzed) {
    cat(sprintf("  • %d voxels excluded (%.1f%%)\n",
                x$metadata$n_voxels_input - x$metadata$n_voxels_analyzed,
                100 * (1 - x$metadata$n_voxels_analyzed/x$metadata$n_voxels_input)))
  }
  
  # Experimental design
  cat("\nDesign:\n")
  cat(sprintf("  • %d conditions\n", x$metadata$n_conditions))
  if (x$metadata$n_trials > 0) {
    cat(sprintf("  • %d trials estimated\n", x$metadata$n_trials))
  }
  
  # HRF information
  cat("\nHRF Estimation:\n")
  cat(sprintf("  • Library: %d HRF variants\n", x$metadata$n_hrfs_library))
  cat(sprintf("  • Manifold: %d dimensions (%s method)\n", 
              x$metadata$manifold_dim, x$metadata$manifold_method))
  cat(sprintf("  • Preset: %s\n", x$metadata$preset_used))
  
  # Quality summary
  if (!is.null(x$qc_metrics)) {
    cat("\nQuality:\n")
    cat(sprintf("  • Mean R²: %.3f\n", x$qc_metrics$mean_r_squared))
    cat(sprintf("  • HRF peak time: %.1f ± %.1f s\n",
                mean(x$qc_metrics$hrf_stats$peak_time, na.rm = TRUE),
                sd(x$qc_metrics$hrf_stats$peak_time, na.rm = TRUE)))
  }
  
  # Runtime
  cat(sprintf("\nProcessing time: %.1f seconds\n", x$metadata$runtime_seconds))
  
  # Hints
  cat("\nUse summary() for detailed statistics\n")
  cat("Use plot() for diagnostic plots\n")
  
  invisible(x)
}


#' Summary method for mhrf_result
#'
#' @param object An mhrf_result object
#' @param ... Additional arguments (ignored)
#' @export
summary.mhrf_result <- function(object, ...) {
  
  # Create summary structure
  summary_obj <- structure(
    list(
      data_dims = c(timepoints = object$metadata$n_timepoints,
                    voxels = object$metadata$n_voxels_analyzed),
      
      design_info = list(
        n_conditions = object$metadata$n_conditions,
        conditions = rownames(object$amplitudes),
        n_trials = object$metadata$n_trials
      ),
      
      hrf_summary = if (!is.null(object$qc_metrics$hrf_stats)) {
        data.frame(
          metric = c("Peak Time (s)", "FWHM (s)", "Time to Peak (s)"),
          mean = c(
            mean(object$qc_metrics$hrf_stats$peak_time, na.rm = TRUE),
            mean(object$qc_metrics$hrf_stats$fwhm, na.rm = TRUE),
            mean(object$qc_metrics$hrf_stats$time_to_peak, na.rm = TRUE)
          ),
          sd = c(
            sd(object$qc_metrics$hrf_stats$peak_time, na.rm = TRUE),
            sd(object$qc_metrics$hrf_stats$fwhm, na.rm = TRUE),
            sd(object$qc_metrics$hrf_stats$time_to_peak, na.rm = TRUE)
          ),
          min = c(
            min(object$qc_metrics$hrf_stats$peak_time, na.rm = TRUE),
            min(object$qc_metrics$hrf_stats$fwhm, na.rm = TRUE),
            min(object$qc_metrics$hrf_stats$time_to_peak, na.rm = TRUE)
          ),
          max = c(
            max(object$qc_metrics$hrf_stats$peak_time, na.rm = TRUE),
            max(object$qc_metrics$hrf_stats$fwhm, na.rm = TRUE),
            max(object$qc_metrics$hrf_stats$time_to_peak, na.rm = TRUE)
          )
        )
      } else NULL,
      
      amplitude_summary = if (!is.null(object$amplitudes) && nrow(object$amplitudes) > 0) {
        # Ensure we have condition names
        if (is.null(rownames(object$amplitudes))) {
          rownames(object$amplitudes) <- paste0("Condition_", 1:nrow(object$amplitudes))
        }
        data.frame(
          condition = rownames(object$amplitudes),
          mean = rowMeans(object$amplitudes),
          sd = apply(object$amplitudes, 1, sd),
          min = apply(object$amplitudes, 1, min),
          max = apply(object$amplitudes, 1, max)
        )
      } else {
        data.frame(
          condition = character(0),
          mean = numeric(0),
          sd = numeric(0),
          min = numeric(0),
          max = numeric(0)
        )
      },
      
      quality_metrics = object$qc_metrics,
      
      parameters = object$metadata$parameters,
      
      object = object
    ),
    class = c("summary.mhrf_result", "list")
  )
  
  return(summary_obj)
}


#' Print method for summary.mhrf_result
#'
#' @param x A summary.mhrf_result object
#' @param ... Additional arguments (ignored)
#' @export
print.summary.mhrf_result <- function(x, ...) {
  cat("\nM-HRF-LSS Analysis Summary\n")
  cat("==========================\n")
  
  # Data dimensions
  cat("\nData Dimensions:\n")
  cat(sprintf("  Timepoints: %d\n", x$data_dims["timepoints"]))
  cat(sprintf("  Voxels analyzed: %d\n", x$data_dims["voxels"]))
  
  # Design
  cat("\nExperimental Design:\n")
  cat(sprintf("  Conditions: %d\n", x$design_info$n_conditions))
  if (x$design_info$n_trials > 0) {
    cat(sprintf("  Trials: %d\n", x$design_info$n_trials))
  }
  
  # HRF characteristics
  if (!is.null(x$hrf_summary)) {
    cat("\nHRF Characteristics:\n")
    print(x$hrf_summary, row.names = FALSE)
  }
  
  # Amplitude summary
  if (nrow(x$amplitude_summary) > 0) {
    cat("\nCondition Amplitudes:\n")
    print(x$amplitude_summary, row.names = FALSE)
  } else {
    cat("\nCondition Amplitudes: None estimated\n")
  }
  
  # Quality metrics
  if (!is.null(x$quality_metrics)) {
    cat("\nQuality Metrics:\n")
    cat(sprintf("  Mean R²: %.3f\n", x$quality_metrics$mean_r_squared))
    cat(sprintf("  Negative amplitudes: %.1f%%\n", 
                x$quality_metrics$percent_negative_amp))
  }
  
  # Key parameters
  cat("\nKey Parameters:\n")
  cat(sprintf("  Preset: %s\n", x$object$metadata$preset_used))
  cat(sprintf("  λ_gamma: %.4f\n", x$parameters$lambda_gamma))
  cat(sprintf("  λ_spatial: %.4f\n", x$parameters$lambda_spatial_smooth))
  cat(sprintf("  Manifold dim: %d\n", x$object$metadata$manifold_dim))
  
  invisible(x)
}


#' Plot method for mhrf_result
#'
#' @param x An mhrf_result object
#' @param type Type of plot: "hrf", "amplitudes", "manifold", "diagnostic"
#' @param voxels Which voxels to plot (for "hrf" type)
#' @param conditions Which conditions to plot (for "amplitudes" type)
#' @param ... Additional graphical parameters
#' @export
plot.mhrf_result <- function(x, type = "diagnostic", voxels = NULL, 
                            conditions = NULL, ...) {
  
  type <- match.arg(type, c("diagnostic", "hrf", "amplitudes", "manifold"))
  
  if (type == "diagnostic") {
    # Create 2x2 diagnostic plot
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    
    # 1. HRF shapes (sample)
    .plot_hrf_shapes(x, n_sample = 20)
    
    # 2. Amplitude distribution
    .plot_amplitude_distribution(x)
    
    # 3. Manifold coordinates
    .plot_manifold_coords(x)
    
    # 4. Quality metrics
    .plot_quality_metrics(x)
    
  } else if (type == "hrf") {
    .plot_hrf_shapes(x, voxels = voxels, ...)
    
  } else if (type == "amplitudes") {
    .plot_amplitude_distribution(x, conditions = conditions, ...)
    
  } else if (type == "manifold") {
    .plot_manifold_coords(x, ...)
  }
  
  invisible(x)
}


# Plotting helper functions ---------------------------------------------------

#' Plot HRF shapes
#' @keywords internal
.plot_hrf_shapes <- function(x, voxels = NULL, n_sample = 20, ...) {
  
  p <- nrow(x$hrf_shapes)
  TR <- x$metadata$parameters$TR
  time_vec <- seq(0, by = TR, length.out = p)
  
  if (is.null(voxels)) {
    # Sample voxels
    n_voxels <- ncol(x$hrf_shapes)
    if (n_voxels > n_sample) {
      voxels <- sample(n_voxels, n_sample)
    } else {
      voxels <- 1:n_voxels
    }
  }
  
  # Plot
  matplot(time_vec, x$hrf_shapes[, voxels], 
          type = "l", lty = 1, col = rgb(0, 0, 0, 0.3),
          xlab = "Time (s)", ylab = "HRF Amplitude",
          main = sprintf("HRF Shapes (n = %d)", length(voxels)))
  
  # Add mean HRF
  mean_hrf <- rowMeans(x$hrf_shapes)
  lines(time_vec, mean_hrf, col = "red", lwd = 3)
  
  # Add canonical for reference if available
  if (requireNamespace("fmrireg", quietly = TRUE)) {
    canonical <- fmrireg::HRF_SPMG1(time_vec)
    canonical <- canonical / sum(abs(canonical))
    lines(time_vec, canonical, col = "blue", lwd = 2, lty = 2)
    legend("topright", c("Estimated", "Mean", "Canonical"), 
           col = c("gray", "red", "blue"), 
           lty = c(1, 1, 2), lwd = c(1, 3, 2))
  } else {
    legend("topright", c("Estimated", "Mean"), 
           col = c("gray", "red"), lty = 1, lwd = c(1, 3))
  }
}


#' Plot amplitude distribution
#' @keywords internal
.plot_amplitude_distribution <- function(x, conditions = NULL, ...) {
  
  if (is.null(conditions)) {
    conditions <- 1:nrow(x$amplitudes)
  }
  
  amplitudes_subset <- x$amplitudes[conditions, , drop = FALSE]
  
  boxplot(t(amplitudes_subset), 
          names = rownames(amplitudes_subset),
          main = "Condition Amplitudes",
          ylab = "Amplitude",
          xlab = "Condition",
          col = rainbow(length(conditions), alpha = 0.5))
  
  abline(h = 0, lty = 2, col = "gray")
}


#' Plot manifold coordinates
#' @keywords internal
.plot_manifold_coords <- function(x, dims = c(1, 2), ...) {
  
  if (nrow(x$manifold_coords) < 2) {
    plot(1, type = "n", main = "Manifold Coordinates",
         xlab = "", ylab = "")
    text(1, 1, "Not enough dimensions to plot", cex = 1.2)
    return()
  }
  
  xi1 <- x$manifold_coords[dims[1], ]
  xi2 <- x$manifold_coords[dims[2], ]
  
  # Color by mean amplitude
  mean_amp <- colMeans(x$amplitudes)
  cols <- colorRampPalette(c("blue", "white", "red"))(100)
  col_idx <- cut(mean_amp, breaks = 100, labels = FALSE)
  
  plot(xi1, xi2, 
       col = cols[col_idx], pch = 19, cex = 0.5,
       main = sprintf("Manifold Coordinates (Dims %d-%d)", dims[1], dims[2]),
       xlab = sprintf("ξ%d", dims[1]),
       ylab = sprintf("ξ%d", dims[2]))
}


#' Plot quality metrics  
#' @keywords internal
.plot_quality_metrics <- function(x, ...) {
  
  if (is.null(x$qc_metrics$hrf_stats)) {
    plot(1, type = "n", main = "Quality Metrics",
         xlab = "", ylab = "")
    text(1, 1, "No quality metrics available", cex = 1.2)
    return()
  }
  
  # Create quality summary plot
  peak_times <- x$qc_metrics$hrf_stats$peak_time
  
  hist(peak_times, breaks = 20,
       main = "HRF Peak Time Distribution",
       xlab = "Peak Time (s)",
       ylab = "Count",
       col = "lightblue")
  
  abline(v = mean(peak_times, na.rm = TRUE), col = "red", lwd = 2)
  abline(v = c(4, 6), col = "green", lty = 2)  # Expected range
  
  legend("topright", 
         c("Mean", "Expected Range"),
         col = c("red", "green"),
         lty = c(1, 2), lwd = c(2, 1))
}


#' Extract coefficients from mhrf_result
#'
#' @param object An mhrf_result object
#' @param type Type of coefficients: "amplitudes", "trial_amplitudes", or "hrfs"
#' @param ... Additional arguments (ignored)
#' @export
coef.mhrf_result <- function(object, type = "amplitudes", ...) {
  
  type <- match.arg(type, c("amplitudes", "trial_amplitudes", "hrfs"))
  
  switch(type,
    amplitudes = object$amplitudes,
    trial_amplitudes = object$trial_amplitudes,
    hrfs = object$hrf_shapes
  )
}


#' Extract fitted values from mhrf_result
#'
#' @param object An mhrf_result object  
#' @param ... Additional arguments (ignored)
#' @export
fitted.mhrf_result <- function(object, ...) {
  if (is.null(object$design_matrices)) {
    stop("Design matrices not available in mhrf_result object")
  }

  X_list <- object$design_matrices
  H <- object$hrf_shapes
  B <- object$amplitudes
  n <- nrow(X_list[[1]])
  V <- ncol(H)
  k <- length(X_list)

  fitted_mat <- matrix(0, n, V)
  for (vx in seq_len(V)) {
    y_hat <- rep(0, n)
    for (c in seq_len(k)) {
      y_hat <- y_hat + (X_list[[c]] %*% H[, vx]) * B[c, vx]
    }
    fitted_mat[, vx] <- y_hat
  }

  return(fitted_mat)
}


#' Convert mhrf_result to data frame
#'
#' @param x An mhrf_result object
#' @param what What to extract: "amplitudes", "hrfs", or "summary"
#' @param ... Additional arguments (ignored)
#' @export
as.data.frame.mhrf_result <- function(x, what = "amplitudes", ...) {
  
  what <- match.arg(what, c("amplitudes", "hrfs", "summary"))
  
  if (what == "amplitudes") {
    # Tidy format for amplitudes
    df <- data.frame(
      voxel = rep(1:ncol(x$amplitudes), each = nrow(x$amplitudes)),
      condition = rep(rownames(x$amplitudes), ncol(x$amplitudes)),
      amplitude = as.vector(x$amplitudes)
    )
    
  } else if (what == "hrfs") {
    # Tidy format for HRFs
    p <- nrow(x$hrf_shapes)
    V <- ncol(x$hrf_shapes)
    TR <- x$metadata$parameters$TR
    
    df <- data.frame(
      voxel = rep(1:V, each = p),
      time = rep(seq(0, by = TR, length.out = p), V),
      hrf = as.vector(x$hrf_shapes)
    )
    
  } else {
    # Summary statistics per voxel
    df <- data.frame(
      voxel = 1:ncol(x$amplitudes),
      mean_amplitude = colMeans(x$amplitudes),
      sd_amplitude = apply(x$amplitudes, 2, sd)
    )
    
    # Add HRF stats if available
    if (!is.null(x$qc_metrics$hrf_stats)) {
      df$hrf_peak_time <- x$qc_metrics$hrf_stats$peak_time
      df$hrf_fwhm <- x$qc_metrics$hrf_stats$fwhm
    }
  }
  
  return(df)
}
