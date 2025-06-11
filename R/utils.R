# Utility functions

#' Default value operator
#'
#' Returns left hand side if not NULL, otherwise right hand side
#' @param a value to check
#' @param b default value
#' @return a if not NULL else b
#' @export
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Adjust HRF Vector for Data Bounds
#'
#' Truncates or pads an HRF vector so that its length does not exceed the
#' available number of time points.
#'
#' @param hrf Numeric vector representing the HRF shape.
#' @param max_timepoints Maximum allowed length.
#' @return HRF vector of length \code{max_timepoints}.
#' @export
adjust_hrf_for_bounds <- function(hrf, max_timepoints) {
  if (!is.numeric(hrf)) {
    stop("hrf must be numeric")
  }
  
  if (!is.numeric(max_timepoints) || length(max_timepoints) != 1 || 
      max_timepoints < 0 || max_timepoints != round(max_timepoints)) {
    stop("max_timepoints must be a non-negative integer")
  }

  if (length(hrf) > max_timepoints) {
    warning("HRF truncated to fit within available timepoints")
    hrf[seq_len(max_timepoints)]
  } else if (length(hrf) < max_timepoints) {
    c(hrf, rep(0, max_timepoints - length(hrf)))
  } else {
    hrf
  }
}


#' Select manifold dimensionality based on eigenvalues
#'
#' Determines the number of diffusion map dimensions needed to explain a
#' desired amount of variance. The function automatically discards the
#' trivial first eigenvalue of the Markov matrix and provides diagnostic
#' logging of the cumulative variance explained and gaps in the spectrum.
#'
#' @param eigenvalues Numeric vector of eigenvalues from the Markov matrix
#'   (ordered by magnitude). The first element is assumed to be the trivial
#'   eigenvalue equal to 1.
#' @param min_var Minimum cumulative variance to retain (between 0 and 1).
#' @return A list with elements:
#'   \itemize{
#'     \item \code{m_auto}: Selected dimensionality.
#'     \item \code{cum_var}: Cumulative variance explained by each dimension.
#'     \item \code{gaps}: Differences between successive eigenvalues.
#'   }
#' @keywords internal
select_manifold_dim <- function(eigenvalues, min_var = 0.95) {
  if (!is.numeric(eigenvalues) || length(eigenvalues) == 0) {
    stop("eigenvalues must be a numeric vector")
  }
  if (min_var <= 0 || min_var > 1) {
    stop("min_var must be between 0 and 1")
  }

  if (length(eigenvalues) == 1) {
    message("Only trivial eigenvalue provided; selecting dimension 1")
    return(list(m_auto = 1, cum_var = 1, gaps = numeric(0)))
  }

  # Exclude the trivial eigenvalue
  eig <- eigenvalues[-1]
  
  # Check for negative eigenvalues 
  # Note: Negative eigenvalues are expected for row-stochastic matrices derived from 
  # Gaussian kernels and do not necessarily indicate numerical errors
  negative_idx <- which(eig < -1e-12)
  
  if (length(negative_idx) > 0) {
    # Only warn for significantly negative eigenvalues
    significant_negative <- sum(eig < -1e-6)
    total_negative <- length(negative_idx)
    
    if (significant_negative > 0) {
      warning(sprintf(
        "Found %d negative eigenvalues (min: %.6e). This is expected for row-stochastic matrices. Using abs() to proceed.",
        total_negative, min(eig[negative_idx])
      ))
    }
    
    # Apply absolute value to all negative eigenvalues
    # This is standard practice in diffusion map literature
    eig[negative_idx] <- abs(eig[negative_idx])
  }
  
  total <- sum(eig)

  if (total <= 0) {
    warning("Non-trivial eigenvalues sum to zero; using dimension 1")
    return(list(m_auto = 1, cum_var = rep(0, length(eig)), gaps = diff(eig)))
  }

  var_explained <- eig / total
  cum_var <- cumsum(var_explained)
  m_auto <- which(cum_var >= min_var)[1]
  if (is.na(m_auto)) {
    m_auto <- length(eig)
    warning(sprintf(
      "Could not achieve %.1f%% variance with available dimensions. Using all %d.",
      min_var * 100, m_auto
    ))
  }

  gaps <- diff(eig)
  gap_idx <- if (length(gaps) > 0) which.max(gaps) + 1 else NA

  msg <- sprintf(
    "Cumulative variance by dimension: %s",
    paste(sprintf("%.3f", cum_var), collapse = " ")
  )
  message(msg)
  if (!is.na(gap_idx)) {
    message(sprintf("Largest eigenvalue gap after dimension %d", gap_idx))
  }

  list(m_auto = m_auto, cum_var = cum_var, gaps = gaps)
}

#' Check RAM feasibility for trial precomputation
#'
#' Estimates expected memory usage for storing trial-by-voxel matrices and
#' compares it to a user-provided limit.
#'
#' @param T_trials Number of trials.
#' @param V Number of voxels.
#' @param ram_limit_GB Memory limit in gigabytes.
#'
#' @return Logical indicating whether precomputation is feasible.
#' @keywords internal
check_ram_feasibility <- function(T_trials, V, ram_limit_GB) {
  if (!is.numeric(T_trials) || length(T_trials) != 1 || T_trials <= 0) {
    stop("T_trials must be a positive numeric scalar")
  }
  if (!is.numeric(V) || length(V) != 1 || V <= 0) {
    stop("V must be a positive numeric scalar")
  }
  if (!is.numeric(ram_limit_GB) || length(ram_limit_GB) != 1 || ram_limit_GB <= 0) {
    stop("ram_limit_GB must be a positive numeric scalar")
  }
  
  expected_GB <- (T_trials * V * 8) / 1e9
  feasible <- expected_GB < ram_limit_GB
  if (!feasible) {
    message(sprintf(
      "Estimated memory %.2f GB exceeds limit %.2f GB - using lazy evaluation for trial regressors.",
      expected_GB, ram_limit_GB
    ))
  }
  feasible
}

#' Validate and standardize ridge penalty parameter
#'
#' Ensures a lambda parameter is a non-negative scalar and applies
#' consistent tolerance-based adjustments. Small values below a fixed
#' threshold are coerced to zero with a warning. Unusually large values
#' trigger a warning about potential over-regularization.
#'
#' @param lambda Numeric value provided by the user.
#' @param name Character name of the parameter (for error messages).
#' @return Sanitized lambda value.
#' @keywords internal
.validate_and_standardize_lambda <- function(lambda, param_name) {
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) {
    stop(param_name, " must be a non-negative scalar")
  }
  as.numeric(lambda)
}
