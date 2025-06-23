# Diagnostic Utility Functions

#' Track Convergence Metrics
#'
#' Tracks convergence metrics across iterations for diagnostics
#'
#' @param current_values Current parameter values (vector or matrix)
#' @param previous_values Previous parameter values
#' @param iteration Current iteration number
#' @param metric_name Name of the metric being tracked
#' @param history Optional existing history to append to
#' @return Updated convergence history
#' @export
track_convergence_metrics <- function(current_values, 
                                    previous_values = NULL,
                                    iteration = 1,
                                    metric_name = "parameters",
                                    history = NULL) {
  
  # Initialize history if needed
  if (is.null(history)) {
    history <- list(
      iterations = numeric(),
      relative_change = numeric(),
      absolute_change = numeric(),
      max_change = numeric(),
      metric_name = metric_name,
      converged = FALSE,
      converged_at = NA
    )
  }
  
  history$iterations <- c(history$iterations, iteration)
  
  # Compute changes if previous values provided
  if (!is.null(previous_values)) {
    # Check for finite values
    finite_mask <- is.finite(current_values) & is.finite(previous_values)
    if (!any(finite_mask)) {
      warning("No finite values for convergence check")
      rel_change <- NA
      abs_change <- NA
      max_change <- NA
    } else {
      current_finite <- current_values[finite_mask]
      previous_finite <- previous_values[finite_mask]
      diff_vals <- current_finite - previous_finite
      
      # Relative change (avoid division by zero and handle NaN/Inf)
      denom <- abs(previous_finite)
      denom[denom < 1e-10] <- 1e-10
      rel_change <- sqrt(mean((diff_vals / denom)^2))
      
      # Absolute change
      abs_change <- sqrt(mean(diff_vals^2))
      
      # Maximum change
      max_change <- max(abs(diff_vals))
    }
    
    history$relative_change <- c(history$relative_change, rel_change)
    history$absolute_change <- c(history$absolute_change, abs_change)
    history$max_change <- c(history$max_change, max_change)
  } else {
    # First iteration
    history$relative_change <- c(history$relative_change, NA)
    history$absolute_change <- c(history$absolute_change, NA)
    history$max_change <- c(history$max_change, NA)
  }
  
  return(history)
}


#' Check Convergence Status
#'
#' Checks if convergence criteria are met and provides diagnostics
#'
#' @param history Convergence history from track_convergence_metrics
#' @param rel_tol Relative tolerance for convergence
#' @param abs_tol Absolute tolerance for convergence
#' @param min_iterations Minimum iterations before checking convergence
#' @param patience Number of iterations to wait for improvement
#' @return List with convergence status and diagnostics
#' @export
check_convergence_status <- function(history,
                                   rel_tol = 1e-4,
                                   abs_tol = 1e-6,
                                   min_iterations = 2,
                                   patience = 3) {
  
  n_iter <- length(history$iterations)
  
  # Need minimum iterations
  if (n_iter < min_iterations) {
    return(list(
      converged = FALSE,
      reason = "insufficient_iterations",
      diagnostic = sprintf("Only %d iterations completed (min: %d)", 
                          n_iter, min_iterations)
    ))
  }
  
  # Get recent changes
  recent_rel <- tail(history$relative_change[!is.na(history$relative_change)], patience)
  recent_abs <- tail(history$absolute_change[!is.na(history$absolute_change)], patience)
  
  # Check relative tolerance
  if (length(recent_rel) >= patience && all(recent_rel < rel_tol)) {
    return(list(
      converged = TRUE,
      reason = "relative_tolerance",
      diagnostic = sprintf("Relative change < %.2e for %d iterations", 
                          rel_tol, patience),
      final_change = tail(recent_rel, 1)
    ))
  }
  
  # Check absolute tolerance
  if (length(recent_abs) >= patience && all(recent_abs < abs_tol)) {
    return(list(
      converged = TRUE,
      reason = "absolute_tolerance",
      diagnostic = sprintf("Absolute change < %.2e for %d iterations", 
                          abs_tol, patience),
      final_change = tail(recent_abs, 1)
    ))
  }
  
  # Check for stagnation (suspicious lack of change)
  if (length(recent_rel) >= 2) {
    rel_var <- var(recent_rel)
    if (rel_var < 1e-10) {
      return(list(
        converged = FALSE,
        reason = "stagnation",
        diagnostic = "Changes are suspiciously constant - may be stuck",
        variance = rel_var
      ))
    }
  }
  
  # Not converged yet
  return(list(
    converged = FALSE,
    reason = "in_progress",
    diagnostic = sprintf("Iteration %d: rel_change = %.2e, abs_change = %.2e",
                        n_iter, 
                        tail(recent_rel, 1),
                        tail(recent_abs, 1))
  ))
}


#' Compute Solution Quality Metrics
#'
#' Computes quality metrics for the current solution
#'
#' @param Y_data Data matrix (n x V)
#' @param Y_predicted Predicted data matrix (n x V)
#' @param Xi_matrix Manifold coordinates (m x V)
#' @param Beta_matrix Amplitude matrix (k x V)
#' @param lambda_smooth Smoothing parameter used
#' @return List of quality metrics
#' @export
compute_solution_quality <- function(Y_data, 
                                   Y_predicted,
                                   Xi_matrix = NULL,
                                   Beta_matrix = NULL,
                                   lambda_smooth = NULL) {
  
  metrics <- list()
  
  # Reconstruction error
  residuals <- Y_data - Y_predicted
  metrics$rmse <- sqrt(mean(residuals^2))
  metrics$r_squared <- 1 - sum(residuals^2) / sum((Y_data - mean(Y_data))^2)
  
  # Per-voxel R-squared
  V <- ncol(Y_data)
  metrics$r_squared_voxels <- numeric(V)
  for (v in 1:V) {
    ss_res <- sum(residuals[, v]^2)
    ss_tot <- sum((Y_data[, v] - mean(Y_data[, v]))^2)
    metrics$r_squared_voxels[v] <- 1 - ss_res / ss_tot
  }
  
  # Smoothness metrics if Xi provided
  if (!is.null(Xi_matrix)) {
    # Spatial smoothness (lower is smoother)
    m <- nrow(Xi_matrix)
    xi_smoothness <- numeric(m)
    for (i in 1:m) {
      # Simple measure: variance of spatial gradients
      xi_row <- Xi_matrix[i, ]
      if (length(xi_row) > 1) {
        xi_smoothness[i] <- var(diff(xi_row))
      }
    }
    metrics$xi_smoothness <- mean(xi_smoothness)
    
    # Effective degrees of freedom (complexity)
    # Approximated by ratio of variance explained
    metrics$effective_df <- sum(apply(Xi_matrix, 1, var) > 1e-6)
  }
  
  # Amplitude statistics if Beta provided
  if (!is.null(Beta_matrix)) {
    metrics$beta_mean <- mean(Beta_matrix)
    metrics$beta_sd <- sd(Beta_matrix)
    metrics$beta_sparsity <- mean(abs(Beta_matrix) < 0.1)
  }
  
  # Overall quality score (0-1)
  quality_components <- c(
    min(metrics$r_squared, 1),  # Fit quality
    1 - min(metrics$xi_smoothness / 10, 1),  # Smoothness (normalized)
    1 - metrics$beta_sparsity  # Non-sparsity
  )
  metrics$overall_quality <- mean(quality_components, na.rm = TRUE)
  
  return(metrics)
}


#' Check Design Matrix Rank
#'
#' Checks rank deficiency in design matrices and removes collinear columns
#'
#' @param X_design Design matrix to check
#' @param tol Tolerance for rank determination
#' @param remove_collinear Whether to remove collinear columns
#' @return List with checked/cleaned design matrix and diagnostics
#' @export
check_design_rank <- function(X_design, tol = 1e-10, remove_collinear = TRUE) {
  
  n <- nrow(X_design)
  p <- ncol(X_design)
  
  # Compute SVD for rank check
  svd_X <- svd(X_design)
  rank_X <- sum(svd_X$d > tol * svd_X$d[1])
  
  result <- list(
    original_rank = rank_X,
    full_rank = p,
    is_rank_deficient = rank_X < p,
    condition_number = svd_X$d[1] / svd_X$d[rank_X]
  )
  
  if (result$is_rank_deficient) {
    warning(sprintf("Design matrix is rank deficient: rank %d < %d columns", 
                   rank_X, p))
    
    if (remove_collinear) {
      # Identify which columns to keep
      # Use QR decomposition with pivoting
      qr_X <- qr(X_design, tol = tol)
      keep_cols <- qr_X$pivot[1:rank_X]
      
      result$X_cleaned <- X_design[, keep_cols, drop = FALSE]
      result$removed_columns <- setdiff(1:p, keep_cols)
      result$kept_columns <- keep_cols
      
      message(sprintf("Removed %d collinear columns: %s", 
                     length(result$removed_columns),
                     paste(result$removed_columns, collapse = ", ")))
    }
  } else {
    result$X_cleaned <- X_design
    result$removed_columns <- integer(0)
    result$kept_columns <- 1:p
  }
  
  return(result)
}


#' Handle Zero Voxels
#'
#' Identifies and handles zero-variance and all-zero voxels
#'
#' @param Y_data Data matrix (n x V)
#' @param min_variance Minimum variance threshold
#' @param replace_with How to handle zero voxels ("skip", "noise", or "mean")
#' @return List with cleaned data and zero voxel indices
#' @export
handle_zero_voxels <- function(Y_data, min_variance = 1e-8, replace_with = "skip") {
  
  n <- nrow(Y_data)
  V <- ncol(Y_data)
  
  # Identify problematic voxels
  # Check for non-finite values first (Inf, -Inf, NaN)
  has_nonfinite <- apply(Y_data, 2, function(x) any(!is.finite(x)))
  
  # Calculate variance safely after checking for finite values
  voxel_vars <- apply(Y_data, 2, function(x) {
    finite_vals <- x[is.finite(x)]
    if (length(finite_vals) < 2) {
      return(NA_real_)  # Not enough finite values
    }
    var(finite_vals)
  })
  
  voxel_means <- colMeans(Y_data, na.rm = TRUE)
  
  all_zero <- abs(voxel_means) < .Machine$double.eps & voxel_vars < .Machine$double.eps
  low_variance <- voxel_vars < min_variance
  
  zero_indices <- which(all_zero | low_variance | has_nonfinite)
  n_zero <- sum(all_zero)
  n_low_var <- sum(low_variance & !all_zero)
  n_nonfinite <- sum(has_nonfinite)
  
  # Report findings
  if (length(zero_indices) > 0) {
    message(sprintf("Found %d problematic voxels: %d all-zero, %d low-variance, %d non-finite", 
                   length(zero_indices), n_zero, n_low_var, n_nonfinite))
  }
  
  # Handle based on strategy
  Y_cleaned <- Y_data
  if (length(zero_indices) > 0 && replace_with != "skip") {
    if (replace_with == "noise") {
      # Replace with small noise
      good_vars <- voxel_vars[voxel_vars > min_variance & is.finite(voxel_vars)]
      if (length(good_vars) > 0) {
        noise_sd <- sqrt(median(good_vars) * 0.01)
      } else {
        noise_sd <- 0.1  # Fallback value
      }
      Y_cleaned[, zero_indices] <- matrix(rnorm(n * length(zero_indices), sd = noise_sd),
                                         n, length(zero_indices))
    } else if (replace_with == "mean") {
      # Replace with global mean signal
      if (length(zero_indices) < V) {
        mean_signal <- rowMeans(Y_data[, -zero_indices, drop = FALSE], na.rm = TRUE)
        Y_cleaned[, zero_indices] <- mean_signal
      }
    }
  }
  
  return(list(
    Y_cleaned = Y_cleaned,
    zero_indices = zero_indices,
    all_zero_indices = which(all_zero),
    low_var_indices = which(low_variance & !all_zero),
    n_problematic = length(zero_indices),
    percent_problematic = 100 * length(zero_indices) / V
  ))
}


#' Monitor Condition Numbers
#'
#' Monitors condition numbers throughout the pipeline
#'
#' @param matrix_list Named list of matrices to check
#' @param warn_threshold Threshold for warning (default 1e8)
#' @param error_threshold Threshold for error (default 1e12)
#' @return Data frame with condition monitoring results
#' @export
monitor_condition_numbers <- function(matrix_list, 
                                    warn_threshold = 1e8,
                                    error_threshold = 1e12) {
  
  results <- data.frame(
    matrix_name = character(),
    dimensions = character(),
    condition_number = numeric(),
    status = character(),
    stringsAsFactors = FALSE
  )
  
  for (name in names(matrix_list)) {
    mat <- matrix_list[[name]]
    
    if (is.matrix(mat) || inherits(mat, "Matrix")) {
      # Compute condition number
      kappa_val <- tryCatch({
        if (nrow(mat) == ncol(mat)) {
          # Square matrix - use standard condition number
          kappa(mat)
        } else {
          # Rectangular - use SVD-based condition
          sv <- svd(mat, nu = 0, nv = 0)
          sv$d[1] / sv$d[length(sv$d)]
        }
      }, error = function(e) NA)
      
      # Determine status
      status <- "OK"
      if (is.na(kappa_val)) {
        status <- "ERROR"
      } else if (kappa_val > error_threshold) {
        status <- "CRITICAL"
        warning(sprintf("Matrix '%s' has critical condition number: %.2e", 
                       name, kappa_val))
      } else if (kappa_val > warn_threshold) {
        status <- "WARNING"
        message(sprintf("Matrix '%s' has high condition number: %.2e", 
                       name, kappa_val))
      }
      
      results <- rbind(results, data.frame(
        matrix_name = name,
        dimensions = sprintf("%dx%d", nrow(mat), ncol(mat)),
        condition_number = kappa_val,
        status = status,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Summary message
  n_warning <- sum(results$status == "WARNING")
  n_critical <- sum(results$status == "CRITICAL")
  
  if (n_critical > 0) {
    warning(sprintf("Found %d matrices with critical conditioning!", n_critical))
  }
  if (n_warning > 0) {
    message(sprintf("Found %d matrices with poor conditioning", n_warning))
  }
  
  return(results)
}


#' Create Progress Bar
#'
#' Creates a simple progress bar for long operations
#'
#' @param total Total number of items
#' @param width Width of progress bar
#' @return Progress bar object
#' @export
create_progress_bar <- function(total, width = 50) {
  
  pb <- list(
    total = total,
    current = 0,
    width = width,
    start_time = Sys.time(),
    last_update = Sys.time()
  )
  
  class(pb) <- "simple_progress_bar"
  return(pb)
}


#' Update Progress Bar
#'
#' Updates and displays progress bar
#'
#' @param pb Progress bar object
#' @param increment Increment amount (default 1)
#' @param message Optional message to display
#' @return Updated progress bar (invisibly)
#' @export
update_progress_bar <- function(pb, increment = 1, message = NULL) {
  
  pb$current <- pb$current + increment
  
  # Only update display every 0.1 seconds to avoid spam
  if (difftime(Sys.time(), pb$last_update, units = "secs") < 0.1 && 
      pb$current < pb$total) {
    return(invisible(pb))
  }
  
  pb$last_update <- Sys.time()
  
  # Calculate progress
  progress <- pb$current / pb$total
  filled <- round(progress * pb$width)
  
  # Time estimation
  elapsed <- difftime(Sys.time(), pb$start_time, units = "secs")
  if (pb$current > 0) {
    eta <- elapsed * (pb$total / pb$current - 1)
    eta_str <- sprintf("ETA: %s", format_time(eta))
  } else {
    eta_str <- "ETA: --:--"
  }
  
  # Build progress bar
  bar <- sprintf("[%s%s] %d/%d (%.1f%%) %s",
                paste(rep("=", filled), collapse = ""),
                paste(rep(" ", pb$width - filled), collapse = ""),
                pb$current, pb$total,
                progress * 100,
                eta_str)
  
  # Add message if provided
  if (!is.null(message)) {
    bar <- paste(bar, "-", message)
  }
  
  # Update display
  cat("\r", bar, sep = "")
  
  # New line when complete
  if (pb$current >= pb$total) {
    cat("\n")
    total_time <- difftime(Sys.time(), pb$start_time, units = "secs")
    message(sprintf("Completed in %s", format_time(total_time)))
  }
  
  return(invisible(pb))
}


#' Format Time Helper
#'
#' Formats seconds into human-readable time
#'
#' @param seconds Number of seconds
#' @return Formatted time string
#' @keywords internal
format_time <- function(seconds) {
  if (seconds < 60) {
    return(sprintf("%.0fs", seconds))
  } else if (seconds < 3600) {
    return(sprintf("%.0fm %.0fs", seconds %/% 60, seconds %% 60))
  } else {
    return(sprintf("%.0fh %.0fm", seconds %/% 3600, (seconds %% 3600) %/% 60))
  }
}
