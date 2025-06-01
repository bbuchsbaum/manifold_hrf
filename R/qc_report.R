# QC Report Generation Functions
# Implementation of MHRF-QC-REPORT-01

#' Generate M-HRF-LSS QC Report
#'
#' Creates a comprehensive HTML quality control report for the M-HRF-LSS pipeline
#' results, including manifold diagnostics, HRF statistics, model fit metrics,
#' and QC flags.
#'
#' @param results An mhrf_results object or list containing pipeline outputs with:
#'   \itemize{
#'     \item \code{core_matrices}: List of core pipeline matrices (Y_data, Xi_smoothed, etc.)
#'     \item \code{manifold}: Manifold construction results
#'     \item \code{diagnostics}: Optional diagnostic metrics (R^2, convergence, etc.)
#'   }
#' @param parameters List of pipeline parameters used, containing:
#'   \itemize{
#'     \item \code{manifold}: Manifold construction parameters
#'     \item \code{pipeline}: Pipeline processing parameters
#'   }
#' @param metadata Optional list of processing metadata (timing, versions, etc.)
#' @param output_file Path to save the HTML report (default: "mhrf_qc_report.html")
#' @param output_dir Directory for report output (default: current directory)
#' @param open_report Logical, whether to open the report in browser (default: TRUE)
#' @param clean_intermediate Logical, whether to remove intermediate files (default: TRUE)
#'
#' @return Path to the generated HTML report
#'
#' @details This function generates a comprehensive QC report including:
#'   \itemize{
#'     \item Executive summary with QC status badges
#'     \item Input parameters table
#'     \item Manifold diagnostics (eigenvalue spectrum, variance explained, reconstruction error)
#'     \item HRF diagnostics (manifold coordinate maps, smoothing effects, shape statistics)
#'     \item Model fit diagnostics (voxel-wise R^2, convergence plots)
#'     \item Performance metrics and timing
#'     \item Failure mode checklist based on QC flags
#'   }
#'
#' @examples
#' \dontrun{
#' # Generate report from pipeline results
#' report_path <- generate_qc_report(
#'   results = pipeline_output,
#'   parameters = list(
#'     manifold = manifold_params,
#'     pipeline = pipeline_params
#'   ),
#'   output_file = "subject001_qc.html"
#' )
#' }
#'
#' @export
generate_qc_report <- function(results,
                              parameters,
                              metadata = NULL,
                              output_file = "mhrf_qc_report.html",
                              output_dir = ".",
                              open_report = TRUE,
                              clean_intermediate = TRUE) {
  
  # Validate inputs
  if (!is.list(results)) {
    stop("results must be a list or mhrf_results object")
  }
  
  if (!is.list(parameters)) {
    stop("parameters must be a list")
  }
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Full output path
  output_path <- file.path(output_dir, output_file)
  
  # Find template
  template_path <- system.file("rmd", "mhrf_qc_report.Rmd", 
                              package = "manifoldhrf")
  
  if (!file.exists(template_path)) {
    stop("QC report template not found. Please reinstall the package.")
  }
  
  # Add diagnostics if not present
  if (is.null(results$diagnostics)) {
    results$diagnostics <- compute_qc_diagnostics(results)
  }
  
  # Prepare parameters for rmarkdown
  render_params <- list(
    results = results,
    parameters = parameters,
    metadata = metadata,
    log = results$log
  )
  
  # Render the report
  message("Generating QC report...")
  
  tryCatch({
    rmarkdown::render(
      input = template_path,
      output_file = basename(output_path),
      output_dir = dirname(output_path),
      params = render_params,
      quiet = TRUE,
      clean = clean_intermediate
    )
    
    message(sprintf("QC report generated: %s", output_path))
    
    # Open in browser if requested
    if (open_report && interactive()) {
      utils::browseURL(output_path)
    }
    
    return(output_path)
    
  }, error = function(e) {
    stop(sprintf("Failed to generate QC report: %s", e$message))
  })
}


#' Compute QC Diagnostics
#' 
#' Computes diagnostic metrics for QC report if not already present
#' 
#' @param results Pipeline results object
#' @return List of diagnostic metrics
#' @keywords internal
compute_qc_diagnostics <- function(results) {
  
  diagnostics <- list()
  
  # Calculate R^2 if we have the necessary components
  if (!is.null(results$core_matrices$Y_data) && 
      !is.null(results$core_matrices$Y_predicted)) {
    
    # Voxel-wise R^2
    Y_data <- results$core_matrices$Y_data
    Y_pred <- results$core_matrices$Y_predicted
    
    diagnostics$r2_voxelwise <- sapply(1:ncol(Y_data), function(v) {
      y <- Y_data[, v]
      y_pred <- Y_pred[, v]
      
      ss_tot <- sum((y - mean(y))^2)
      ss_res <- sum((y - y_pred)^2)
      
      if (ss_tot > 0) {
        1 - ss_res / ss_tot
      } else {
        NA
      }
    })
  }
  
  # Calculate reconstruction error if we have manifold components
  if (!is.null(results$manifold$L_library) && 
      !is.null(results$manifold$B_reconstructor) &&
      !is.null(results$manifold$Phi_coords)) {
    
    L <- results$manifold$L_library
    B <- results$manifold$B_reconstructor
    Phi <- results$manifold$Phi_coords
    
    # Calculate error for different dimensions
    max_dim <- ncol(B)
    diagnostics$reconstruction_error <- numeric(max_dim)
    
    for (m in 1:max_dim) {
      L_recon <- B[, 1:m] %*% t(Phi[, 1:m])
      error <- norm(L - L_recon, "F") / norm(L, "F")
      diagnostics$reconstruction_error[m] <- error
    }
  }
  
  # Add convergence tracking if alternating optimization was used
  if (!is.null(results$convergence_history)) {
    diagnostics$convergence <- results$convergence_history
  }
  
  return(diagnostics)
}


#' Create QC Flag Summary
#'
#' Evaluates pipeline results and generates QC flags for potential issues
#'
#' @param results Pipeline results object
#' @param thresholds List of threshold values for QC checks
#' @return List of QC flags with status and messages
#' @export
create_qc_flags <- function(results, 
                           thresholds = list(
                             min_trials = 20,
                             min_r2 = 0.1,
                             max_poor_fit_fraction = 0.3,
                             min_hrf_peak = 2,
                             max_hrf_peak = 10,
                             max_motion_dvars = 1.5
                           )) {
  
  qc_flags <- list()
  
  # Check trial count
  if (!is.null(results$core_matrices$Beta_trial)) {
    n_trials <- nrow(results$core_matrices$Beta_trial)
    if (n_trials < thresholds$min_trials) {
      qc_flags$low_trial_count <- list(
        status = "warning",
        message = sprintf("Low trial count: %d trials (recommended: ≥%d)", 
                         n_trials, thresholds$min_trials),
        severity = 2
      )
    }
  }
  
  # Check model fit
  if (!is.null(results$diagnostics$r2_voxelwise)) {
    poor_fits <- mean(results$diagnostics$r2_voxelwise < thresholds$min_r2, 
                     na.rm = TRUE)
    
    if (poor_fits > thresholds$max_poor_fit_fraction) {
      qc_flags$poor_fits <- list(
        status = "warning",
        message = sprintf("%.1f%% of voxels have R^2 < %.2f", 
                         100 * poor_fits, thresholds$min_r2),
        severity = 2
      )
    }
  }
  
  # Check HRF characteristics
  if (!is.null(results$hrf_stats)) {
    unusual_peaks <- sum(
      results$hrf_stats$peak_time < thresholds$min_hrf_peak |
      results$hrf_stats$peak_time > thresholds$max_hrf_peak,
      na.rm = TRUE
    )
    
    if (unusual_peaks > 0) {
      n_voxels <- length(results$hrf_stats$peak_time)
      qc_flags$unstable_hrf <- list(
        status = "warning",
        message = sprintf("%d voxels (%.1f%%) have unusual HRF peaks (<%ds or >%ds)",
                         unusual_peaks, 100 * unusual_peaks / n_voxels,
                         thresholds$min_hrf_peak, thresholds$max_hrf_peak),
        severity = 3
      )
    }
  }
  
  # Check motion if available
  if (!is.null(results$metadata$motion_dvars)) {
    high_motion_frames <- sum(results$metadata$motion_dvars > thresholds$max_motion_dvars)
    if (high_motion_frames > 0) {
      qc_flags$high_motion <- list(
        status = "warning",
        message = sprintf("%d frames with high motion (DVARS > %.1f)",
                         high_motion_frames, thresholds$max_motion_dvars),
        severity = 2
      )
    }
  }
  
  # Overall QC status
  if (length(qc_flags) == 0) {
    qc_flags$overall <- list(
      status = "pass",
      message = "All QC checks passed",
      severity = 0
    )
  } else {
    max_severity <- max(sapply(qc_flags, function(x) x$severity))
    qc_flags$overall <- list(
      status = ifelse(max_severity >= 3, "fail", "warning"),
      message = sprintf("%d QC issues detected", length(qc_flags) - 1),
      severity = max_severity
    )
  }
  
  class(qc_flags) <- c("mhrf_qc_flags", "list")
  return(qc_flags)
}


#' Print QC Flags
#' @export
print.mhrf_qc_flags <- function(x, ...) {
  cat("M-HRF-LSS QC Flags Summary\n")
  cat("==========================\n\n")
  
  # Overall status
  overall <- x$overall
  status_symbol <- switch(overall$status,
    "pass" = "✓",
    "warning" = "⚠",
    "fail" = "✗"
  )
  
  cat(sprintf("Overall Status: %s %s\n\n", status_symbol, toupper(overall$status)))
  
  # Individual flags
  for (name in names(x)) {
    if (name != "overall") {
      flag <- x[[name]]
      # Convert underscore to space and capitalize each word
      formatted_name <- paste0(toupper(substring(strsplit(name, "_")[[1]], 1, 1)),
                              substring(strsplit(name, "_")[[1]], 2))
      formatted_name <- paste(formatted_name, collapse = " ")
      cat(sprintf("- %s: %s\n", formatted_name, flag$message))
    }
  }
}


#' Extract HRF Statistics
#'
#' Computes statistics about estimated HRF shapes
#'
#' @param H_shapes Matrix of HRF shapes (p x V)
#' @param TR_precision Time resolution of HRF sampling
#' @return Data frame with HRF statistics per voxel
#' @export
extract_hrf_stats <- function(H_shapes, TR_precision = 0.1) {
  
  if (!is.matrix(H_shapes)) {
    stop("H_shapes must be a matrix")
  }
  
  p <- nrow(H_shapes)
  V <- ncol(H_shapes)
  
  # Time grid
  time_points <- seq(0, (p - 1) * TR_precision, by = TR_precision)
  
  # Initialize results
  hrf_stats <- data.frame(
    voxel = 1:V,
    peak_time = NA_real_,
    peak_amplitude = NA_real_,
    time_to_peak = NA_real_,
    fwhm = NA_real_,
    undershoot_ratio = NA_real_
  )
  
  for (v in 1:V) {
    hrf <- H_shapes[, v]
    
    # Skip if HRF is all zeros or NAs
    if (all(is.na(hrf)) || all(hrf == 0)) {
      next
    }
    
    # Find peak
    peak_idx <- which.max(hrf)
    hrf_stats$peak_time[v] <- time_points[peak_idx]
    hrf_stats$peak_amplitude[v] <- hrf[peak_idx]
    hrf_stats$time_to_peak[v] <- time_points[peak_idx]
    
    # Estimate FWHM
    if (max(hrf) > 0) {
      half_max <- max(hrf) / 2
      above_half <- which(hrf >= half_max)
      
      if (length(above_half) >= 2) {
        hrf_stats$fwhm[v] <- time_points[max(above_half)] - time_points[min(above_half)]
      }
    }
    
    # Undershoot ratio (if present)
    if (peak_idx < length(hrf)) {
      post_peak <- hrf[(peak_idx + 1):length(hrf)]
      if (length(post_peak) > 0 && any(post_peak < 0) && !is.na(hrf[peak_idx]) && hrf[peak_idx] != 0) {
        hrf_stats$undershoot_ratio[v] <- abs(min(post_peak)) / hrf[peak_idx]
      }
    }
  }
  
  return(hrf_stats)
}
