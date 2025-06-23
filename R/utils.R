# Utility functions

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
#' @param verbose Logical whether to show diagnostic messages about negative
#'   eigenvalues (default: FALSE).
#' @return A list with elements:
#'   \itemize{
#'     \item \code{m_auto}: Selected dimensionality.
#'     \item \code{cum_var}: Cumulative variance explained by each dimension.
#'     \item \code{gaps}: Differences between successive eigenvalues.
#'   }
#' @keywords internal
select_manifold_dim <- function(eigenvalues, min_var = 0.95, verbose = FALSE) {
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
    # Only log for significantly negative eigenvalues
    significant_negative <- sum(eig < -1e-6)
    total_negative <- length(negative_idx)
    
    if (significant_negative > 0 && verbose) {
      # Log as debug message only if verbose mode is on
      # Negative eigenvalues are normal for row-stochastic matrices from Gaussian kernels
      message(sprintf(
        "Note: Found %d negative eigenvalues (min: %.6e) in diffusion map. This is normal for row-stochastic matrices.",
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
  
  # Convert to numeric to avoid integer overflow
  # Calculate memory in GB directly to avoid overflow
  T_trials <- as.numeric(T_trials)
  V <- as.numeric(V)
  
  # Calculate GB step by step to avoid overflow
  # 8 bytes per double, 1e9 bytes per GB
  expected_GB <- (T_trials / 1e9) * V * 8
  
  # If still too large, try alternate calculation
  if (!is.finite(expected_GB)) {
    expected_GB <- T_trials * (V / 1e9) * 8
  }
  
  # If still not finite, estimate is too large
  if (!is.finite(expected_GB)) {
    message("Memory requirement exceeds computational limits - using lazy evaluation for trial regressors.")
    return(FALSE)
  }
  
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

#' Check Memory Requirements for Distance Matrix
#'
#' Checks if creating a distance matrix would exceed available memory.
#' 
#' @param n_points Number of points in the distance matrix
#' @param element_size Size of each element in bytes (default: 8 for double)
#' @param safety_factor Fraction of available memory to use (default: 0.5)
#' @return NULL (invisibly). Throws error if memory would be exceeded.
#' @keywords internal
check_distance_memory <- function(n_points, element_size = 8, safety_factor = 0.5) {
  # Calculate required memory in GB
  estimated_gb <- (n_points * n_points * element_size) / 1e9
  
  # Get available memory
  available_gb <- .get_available_memory_gb()
  
  # Check if we would exceed memory limits
  if (estimated_gb > available_gb * safety_factor) {
    stop(sprintf(
      "Distance matrix would require %.1f GB but only %.1f GB available (safety factor: %.0f%%)\n",
      estimated_gb, available_gb, safety_factor * 100),
      "Consider using sparse methods or processing in chunks.",
      call. = FALSE)
  }
  
  invisible(NULL)
}

#' Get Available Memory in GB
#' 
#' Platform-agnostic function to get available memory.
#' 
#' @return Available memory in GB
#' @keywords internal
.get_available_memory_gb <- function() {
  # Final fallback if all methods fail
  fallback_gb <- 8
  
  mem_gb <- tryCatch({
    sysname <- Sys.info()["sysname"]
    
    # Use OS-specific logic
    switch(sysname,
      "Linux" = {
        # Primary: Use MemAvailable from /proc/meminfo (in KB)
        meminfo <- readLines("/proc/meminfo")
        available_line <- grep("^MemAvailable:", meminfo, value = TRUE)
        if (length(available_line) > 0) {
          kb <- as.numeric(gsub("[^0-9]", "", available_line))
          return(kb / 1024^2)  # KB to GB
        }
        
        # Fallback: Use MemTotal if MemAvailable isn't present
        total_line <- grep("^MemTotal:", meminfo, value = TRUE)
        if (length(total_line) > 0) {
          warning("Could not determine available memory on Linux, using total memory", call. = FALSE)
          kb <- as.numeric(gsub("[^0-9]", "", total_line))
          return(kb / 1024^2)
        }
        
        stop("Could not parse /proc/meminfo")
      },
      
      "Darwin" = {
        # Primary: Calculate available memory from vm_stat
        tryCatch({
          pagesize <- as.numeric(system("sysctl -n hw.pagesize", intern = TRUE))
          vm_stat <- system("vm_stat", intern = TRUE)
          free_pages <- as.numeric(gsub("[^0-9]", "", grep("Pages free", vm_stat, value = TRUE)))
          inactive_pages <- as.numeric(gsub("[^0-9]", "", grep("Pages inactive", vm_stat, value = TRUE)))
          # Free + inactive pages gives good estimate of available memory
          bytes <- (free_pages + inactive_pages) * pagesize
          return(bytes / 1024^3)  # Bytes to GB
        }, error = function(e) {
          # Fallback: Get total memory from sysctl
          warning("Could not determine available memory on macOS, using total memory", call. = FALSE)
          bytes <- as.numeric(system("sysctl -n hw.memsize", intern = TRUE))
          return(bytes / 1024^3)
        })
      },
      
      "Windows" = {
        # Primary: Use wmic to get FreePhysicalMemory (in KB)
        tryCatch({
          mem_out <- system("wmic OS get FreePhysicalMemory /Value", intern = TRUE)
          val_line <- grep("=", mem_out, value = TRUE)
          if (length(val_line) > 0) {
            kb <- as.numeric(trimws(strsplit(val_line, "=")[[1]][2]))
            return(kb / 1024^2)
          }
          stop("Could not parse wmic output")
        }, error = function(e) {
          # Fallback: Get TotalPhysicalMemory (in Bytes)
          warning("Could not determine free memory on Windows, using total memory", call. = FALSE)
          tryCatch({
            mem_out <- system("wmic ComputerSystem get TotalPhysicalMemory /Value", intern = TRUE)
            val_line <- grep("=", mem_out, value = TRUE)
            if (length(val_line) > 0) {
              bytes <- as.numeric(trimws(strsplit(val_line, "=")[[1]][2]))
              return(bytes / 1024^3)
            }
          }, error = function(e2) {
            # Last resort for Windows - try memory.limit()
            if (exists("memory.limit")) {
              return(memory.limit() / 1024)  # MB to GB
            }
          })
          stop("All Windows memory detection methods failed")
        })
      },
      
      # Unsupported OS
      stop(paste("Unsupported operating system:", sysname))
    )
  }, error = function(e) {
    # Final catch-all
    warning(paste0("Could not determine system memory, assuming ", fallback_gb, "GB. Error: ", e$message), 
            call. = FALSE)
    return(fallback_gb)
  })
  
  # Ensure positive value
  max(mem_gb, 1, na.rm = TRUE)
}
