# Input Validation and Error Handling for M-HRF-LSS
# Provides comprehensive validation with helpful error messages

#' Validate Y_data Input
#'
#' @param Y_data Input data matrix or neuroimaging object
#' @return Validated data information
#' @keywords internal
.validate_Y_data <- function(Y_data) {
  
  # Check basic type
  if (is.null(Y_data)) {
    stop("Y_data cannot be NULL", call. = FALSE)
  }
  
  if (!is.matrix(Y_data) && !inherits(Y_data, c("NeuroVec", "NeuroVol"))) {
    stop(
      "Y_data must be a matrix or NeuroVec/NeuroVol object\n",
      "  Received: ", class(Y_data)[1], "\n",
      "  Hint: Use as.matrix() to convert data frames to matrices",
      call. = FALSE
    )
  }
  
  # Convert to matrix for validation
  if (inherits(Y_data, c("NeuroVec", "NeuroVol"))) {
    if (!requireNamespace("neuroim2", quietly = TRUE)) {
      stop(
        "Package 'neuroim2' required for neuroimaging data\n",
        "  Install with: install.packages('neuroim2')",
        call. = FALSE
      )
    }
    Y_matrix <- as.matrix(Y_data)
  } else {
    Y_matrix <- Y_data
  }
  
  # Check dimensions
  if (nrow(Y_matrix) < 10) {
    stop(
      "Y_data has too few timepoints: ", nrow(Y_matrix), "\n",
      "  Minimum required: 10 timepoints\n",
      "  Hint: Each row should be a timepoint, each column a voxel",
      call. = FALSE
    )
  }
  
  if (ncol(Y_matrix) < 1) {
    stop(
      "Y_data has no voxels\n",
      "  Hint: Each column should represent one voxel",
      call. = FALSE
    )
  }
  
  # Check for problematic values
  if (any(!is.finite(Y_matrix))) {
    n_bad <- sum(!is.finite(Y_matrix))
    total <- length(Y_matrix)
    pct_bad <- 100 * n_bad / total
    
    if (pct_bad > 50) {
      stop(
        "Y_data contains too many non-finite values: ", round(pct_bad, 1), "%\n",
        "  Found ", n_bad, " non-finite values out of ", total, " total\n",
        "  Hint: Check for NaN, Inf, or NA values in your data",
        call. = FALSE
      )
    } else if (pct_bad > 10) {
      warning(
        "Y_data contains ", round(pct_bad, 1), "% non-finite values\n",
        "  This may affect analysis quality"
      )
    }
  }
  
  # Check variance
  voxel_vars <- apply(Y_matrix, 2, var, na.rm = TRUE)
  zero_var_voxels <- sum(voxel_vars < .Machine$double.eps, na.rm = TRUE)
  
  if (zero_var_voxels == ncol(Y_matrix)) {
    stop(
      "All voxels have zero variance\n",
      "  Hint: Check if Y_data contains actual signal",
      call. = FALSE
    )
  }
  
  if (zero_var_voxels > ncol(Y_matrix) * 0.5) {
    warning(
      "Many voxels have zero variance: ", zero_var_voxels, " out of ", ncol(Y_matrix), "\n",
      "  Consider using a brain mask to exclude non-brain voxels"
    )
  }
  
  return(list(
    n_timepoints = nrow(Y_matrix),
    n_voxels = ncol(Y_matrix),
    has_signal = zero_var_voxels < ncol(Y_matrix) * 0.9
  ))
}


#' Validate Events Data Frame
#'
#' @param events Events data frame
#' @param n_timepoints Number of timepoints in Y_data
#' @param TR Repetition time
#' @return Validated events information
#' @keywords internal
.validate_events <- function(events, n_timepoints, TR) {
  
  # Check basic type
  if (is.null(events)) {
    stop("events cannot be NULL", call. = FALSE)
  }
  
  if (!is.data.frame(events)) {
    stop(
      "events must be a data frame\n",
      "  Received: ", class(events)[1], "\n",
      "  Hint: Use data.frame() to create events table",
      call. = FALSE
    )
  }
  
  # Check for empty data frame
  if (nrow(events) == 0) {
    stop(
      "events data frame is empty\n",
      "  Hint: Make sure your events data frame has at least one row",
      call. = FALSE
    )
  }
  
  # Check required columns
  required_cols <- c("condition", "onset")
  missing_cols <- setdiff(required_cols, names(events))
  
  if (length(missing_cols) > 0) {
    stop(
      "events data frame missing required columns: ", paste(missing_cols, collapse = ", "), "\n",
      "  Required columns: ", paste(required_cols, collapse = ", "), "\n",
      "  Found columns: ", paste(names(events), collapse = ", "), "\n",
      "  Hint: Ensure your events data frame has 'condition' and 'onset' columns",
      call. = FALSE
    )
  }
  
  # Add duration if missing
  if (!"duration" %in% names(events)) {
    events$duration <- 0
    message("Added default duration = 0 for all events")
  }
  
  # Validate onset times
  if (any(!is.numeric(events$onset))) {
    stop(
      "onset column must be numeric (time in seconds)\n",
      "  Hint: Convert onset times to numeric values",
      call. = FALSE
    )
  }
  
  if (any(events$onset < 0)) {
    stop(
      "onset times cannot be negative\n",
      "  Found negative onsets: ", paste(which(events$onset < 0), collapse = ", "),
      call. = FALSE
    )
  }
  
  max_time <- (n_timepoints - 1) * TR
  late_events <- events$onset > max_time
  
  if (any(late_events)) {
    n_late <- sum(late_events)
    stop(
      "Some events occur after data ends\n",
      "  Data ends at: ", max_time, " seconds\n",
      "  Late events: ", n_late, " (rows: ", paste(which(late_events), collapse = ", "), ")\n",
      "  Hint: Check TR and make sure onset times match your data length",
      call. = FALSE
    )
  }
  
  # Validate duration
  if (any(!is.numeric(events$duration))) {
    stop(
      "duration column must be numeric (duration in seconds)\n",
      "  Hint: Convert durations to numeric values",
      call. = FALSE
    )
  }
  
  if (any(events$duration < 0)) {
    stop(
      "duration cannot be negative\n",
      "  Found negative durations: ", paste(which(events$duration < 0), collapse = ", "),
      call. = FALSE
    )
  }
  
  # Validate conditions
  if (any(is.na(events$condition))) {
    stop(
      "condition column cannot contain NA values\n",
      "  Found NA conditions in rows: ", paste(which(is.na(events$condition)), collapse = ", "),
      call. = FALSE
    )
  }
  
  conditions <- unique(events$condition)
  n_conditions <- length(conditions)
  
  if (n_conditions > 20) {
    warning(
      "Many conditions detected: ", n_conditions, "\n",
      "  This may lead to overfitting. Consider grouping similar conditions."
    )
  }
  
  # Check event spacing
  events_sorted <- events[order(events$onset), ]
  if (nrow(events_sorted) > 1) {
    gaps <- diff(events_sorted$onset)
    min_gap <- min(gaps)
    
    if (min_gap < TR) {
      warning(
        "Some events are very close together (min gap: ", round(min_gap, 2), " s)\n",
        "  This may cause collinearity issues"
      )
    }
  }
  
  return(list(
    n_events = nrow(events),
    n_conditions = n_conditions,
    conditions = conditions,
    has_trials = "trial_type" %in% names(events) || "trial_id" %in% names(events),
    events_validated = events
  ))
}


#' Validate Analysis Parameters
#'
#' @param TR Repetition time
#' @param preset Parameter preset name
#' @param n_voxels Number of voxels
#' @param user_params Additional user parameters
#' @return Validated parameters
#' @keywords internal
.validate_parameters <- function(TR, preset, n_voxels, user_params = list()) {
  
  # Validate TR
  if (!is.numeric(TR) || length(TR) != 1) {
    stop(
      "TR must be a single numeric value\n",
      "  Received: ", class(TR)[1], " of length ", length(TR), "\n",
      "  Hint: TR should be repetition time in seconds (e.g., TR = 2)",
      call. = FALSE
    )
  }
  
  if (TR <= 0 || TR > 10) {
    stop(
      "TR value seems unrealistic: ", TR, " seconds\n",
      "  Expected range: 0.1 to 10 seconds\n",
      "  Hint: Make sure TR is in seconds, not milliseconds",
      call. = FALSE
    )
  }
  
  # Validate preset
  valid_presets <- c("conservative", "balanced", "aggressive", "fast", "quality", "robust")
  if (!preset %in% valid_presets) {
    stop(
      "Invalid preset: '", preset, "'\n",
      "  Valid presets: ", paste(valid_presets, collapse = ", "), "\n",
      "  Hint: Use 'balanced' for most analyses",
      call. = FALSE
    )
  }
  
  # Check for large datasets
  if (n_voxels > 100000) {
    warning(
      "Large dataset detected: ", n_voxels, " voxels\n",
      "  Consider using chunked processing or a brain mask\n",
      "  Hint: Use preset = 'fast' for faster processing"
    )
  }
  
  # Validate common user parameters
  if ("lambda_gamma" %in% names(user_params)) {
    lambda_gamma <- user_params$lambda_gamma
    if (!is.numeric(lambda_gamma) || lambda_gamma < 0) {
      stop(
        "lambda_gamma must be a non-negative number\n",
        "  Received: ", lambda_gamma, "\n",
        "  Hint: Use values between 0.001 and 0.1",
        call. = FALSE
      )
    }
    if (lambda_gamma > 1) {
      warning("lambda_gamma is very large: ", lambda_gamma, ". This may cause over-regularization.")
    }
  }
  
  if ("m_manifold_dim_target" %in% names(user_params)) {
    m_dim <- user_params$m_manifold_dim_target
    if (!is.numeric(m_dim) || m_dim < 1 || m_dim != round(m_dim)) {
      stop(
        "m_manifold_dim_target must be a positive integer\n",
        "  Received: ", m_dim, "\n",
        "  Hint: Use values between 3 and 10",
        call. = FALSE
      )
    }
    if (m_dim > 15) {
      warning("m_manifold_dim_target is very large: ", m_dim, ". This may cause overfitting.")
    }
  }
  
  return(TRUE)
}


#' Validate Voxel Mask
#'
#' @param voxel_mask Voxel mask (logical or numeric)
#' @param n_voxels Total number of voxels
#' @return Validated mask information
#' @keywords internal
.validate_voxel_mask <- function(voxel_mask, n_voxels) {
  
  if (is.null(voxel_mask)) {
    return(list(
      use_mask = FALSE,
      n_voxels_kept = n_voxels
    ))
  }
  
  # Check type
  if (!is.logical(voxel_mask) && !is.numeric(voxel_mask)) {
    stop(
      "voxel_mask must be logical or numeric\n",
      "  Received: ", class(voxel_mask)[1], "\n",
      "  Hint: Use TRUE/FALSE values or 0/1 values",
      call. = FALSE
    )
  }
  
  # Check length
  if (length(voxel_mask) != n_voxels) {
    stop(
      "voxel_mask length doesn't match number of voxels\n",
      "  Mask length: ", length(voxel_mask), "\n",
      "  Number of voxels: ", n_voxels, "\n",
      "  Hint: Mask should have one element per voxel",
      call. = FALSE
    )
  }
  
  # Convert to logical if numeric
  if (is.numeric(voxel_mask)) {
    if (any(!voxel_mask %in% c(0, 1))) {
      warning("Numeric voxel_mask contains values other than 0 and 1. Using > 0 as criterion.")
    }
    voxel_mask <- voxel_mask > 0
  }
  
  n_kept <- sum(voxel_mask)
  
  if (n_kept == 0) {
    stop(
      "voxel_mask excludes all voxels\n",
      "  Hint: Make sure your mask has some TRUE values",
      call. = FALSE
    )
  }
  
  if (n_kept < 10) {
    warning(
      "Very few voxels selected: ", n_kept, "\n",
      "  This may not provide sufficient data for reliable analysis"
    )
  }
  
  pct_kept <- 100 * n_kept / n_voxels
  message(sprintf("Voxel mask: keeping %d/%d voxels (%.1f%%)", n_kept, n_voxels, pct_kept))
  
  return(list(
    use_mask = TRUE,
    mask = voxel_mask,
    n_voxels_kept = n_kept,
    percent_kept = pct_kept
  ))
}


#' Check System Requirements
#'
#' @param n_voxels Number of voxels
#' @param n_timepoints Number of timepoints
#' @param preset Analysis preset
#' @return System check results
#' @keywords internal
.check_system_requirements <- function(n_voxels, n_timepoints, preset) {
  
  # Estimate memory requirements
  data_size_gb <- n_voxels * n_timepoints * 8 / (1024^3)  # 8 bytes per double
  
  # Get available memory (rough estimate)
  available_gb <- tryCatch({
    if (Sys.info()["sysname"] == "Linux") {
      # Try to read from /proc/meminfo
      meminfo <- readLines("/proc/meminfo", n = 10)
      available_line <- grep("MemAvailable", meminfo, value = TRUE)
      if (length(available_line) > 0) {
        available_kb <- as.numeric(gsub("\\D", "", available_line))
        available_kb / (1024^2)
      } else {
        NA
      }
    } else {
      NA  # Can't easily determine on other systems
    }
  }, error = function(e) NA)
  
  warnings <- character(0)
  
  # Memory warnings
  if (data_size_gb > 2) {
    warnings <- c(warnings, 
      sprintf("Large dataset: %.1f GB of data", data_size_gb))
    
    if (!is.na(available_gb) && data_size_gb > available_gb * 0.5) {
      warnings <- c(warnings,
        "Dataset may exceed available memory. Consider using chunked processing.")
    }
  }
  
  # Performance warnings
  if (n_voxels > 50000 && preset %in% c("quality", "aggressive")) {
    warnings <- c(warnings,
      "Large dataset with intensive preset may be very slow. Consider 'fast' or 'balanced' preset.")
  }
  
  if (length(warnings) > 0) {
    for (w in warnings) {
      warning(w, call. = FALSE, immediate. = TRUE)
    }
  }
  
  return(list(
    data_size_gb = data_size_gb,
    warnings = warnings,
    estimated_time_minutes = if (n_voxels > 10000) n_voxels / 1000 else 1
  ))
}