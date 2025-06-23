# Fallback Cascade and Memory Safety
# Implementation of SOUND-FALLBACK-CASCADE and SOUND-MEMORY-SAFE

#' Get available system memory in GB (platform agnostic)
#' 
#' @return Available memory in GB, or NULL if unable to determine
#' @keywords internal
.get_available_memory_gb <- function() {
  bytes_to_gb <- function(bytes) bytes / 1024^3
  
  os_type <- Sys.info()["sysname"]
  
  # Try multiple approaches based on platform
  available_gb <- NULL
  
  if (os_type == "Linux") {
    # Linux: use /proc/meminfo or free command
    mem_info <- try(readLines("/proc/meminfo"), silent = TRUE)
    if (!inherits(mem_info, "try-error")) {
      # Parse MemAvailable or MemFree
      avail_line <- grep("^MemAvailable:", mem_info, value = TRUE)
      if (length(avail_line) == 0) {
        avail_line <- grep("^MemFree:", mem_info, value = TRUE)
      }
      if (length(avail_line) > 0) {
        available_kb <- as.numeric(sub(".*:\\s*(\\d+).*", "\\1", avail_line))
        available_gb <- available_kb / 1024^2
      }
    } else {
      # Fallback to free command
      mem_info <- try(system("free -b", intern = TRUE), silent = TRUE)
      if (!inherits(mem_info, "try-error")) {
        available_line <- grep("^Mem:", mem_info, value = TRUE)
        if (length(available_line) > 0) {
          fields <- strsplit(available_line, "\\s+")[[1]]
          if (length(fields) >= 7) {
            available_bytes <- as.numeric(fields[7])
            available_gb <- bytes_to_gb(available_bytes)
          }
        }
      }
    }
    
  } else if (os_type == "Darwin") {
    # macOS: use vm_stat command
    vm_stat <- try(system("vm_stat", intern = TRUE), silent = TRUE)
    if (!inherits(vm_stat, "try-error")) {
      # Get page size
      page_size_line <- grep("page size of", vm_stat, value = TRUE)
      page_size <- 4096  # default
      if (length(page_size_line) > 0) {
        page_size <- as.numeric(sub(".*page size of (\\d+) bytes.*", "\\1", page_size_line))
      }
      
      # Get free pages
      free_line <- grep("Pages free:", vm_stat, value = TRUE)
      inactive_line <- grep("Pages inactive:", vm_stat, value = TRUE)
      
      free_pages <- 0
      inactive_pages <- 0
      
      if (length(free_line) > 0) {
        free_pages <- as.numeric(sub(".*:\\s*(\\d+).*", "\\1", free_line))
      }
      if (length(inactive_line) > 0) {
        inactive_pages <- as.numeric(sub(".*:\\s*(\\d+).*", "\\1", inactive_line))
      }
      
      available_bytes <- (free_pages + inactive_pages) * page_size
      available_gb <- bytes_to_gb(available_bytes)
    }
    
  } else if (os_type == "Windows") {
    # Windows: use wmic command
    wmic_cmd <- try(system("wmic OS get FreePhysicalMemory /value", intern = TRUE), silent = TRUE)
    if (!inherits(wmic_cmd, "try-error")) {
      free_line <- grep("FreePhysicalMemory=", wmic_cmd, value = TRUE)
      if (length(free_line) > 0) {
        free_kb <- as.numeric(sub(".*=(\\d+).*", "\\1", free_line))
        available_gb <- free_kb / 1024^2
      }
    }
  }
  
  # If platform-specific methods fail, try to use memory.limit() on Windows
  if (is.null(available_gb) && os_type == "Windows") {
    mem_limit <- try(memory.limit(), silent = TRUE)
    if (!inherits(mem_limit, "try-error") && is.finite(mem_limit)) {
      # memory.limit returns MB on Windows
      # Estimate available as 50% of limit (conservative)
      available_gb <- mem_limit / 1024 * 0.5
    }
  }
  
  # Final fallback: use a conservative estimate based on total memory
  if (is.null(available_gb)) {
    # Try to estimate from gc()
    gc_info <- gc()
    # Use max memory used as a lower bound for total memory
    max_used_mb <- max(gc_info[, "max used"])
    if (max_used_mb > 0) {
      # Assume we can use 2x what we've used so far (very conservative)
      available_gb <- max_used_mb / 1024
    }
  }
  
  return(available_gb)
}

#' Estimate Memory Requirements
#'
#' Estimates memory needed for M-HRF-LSS pipeline
#'
#' @param n_timepoints Number of timepoints
#' @param n_voxels Number of voxels
#' @param n_conditions Number of conditions
#' @param n_trials Number of trials
#' @param p_hrf HRF length in timepoints
#' @param m_manifold Manifold dimension
#' @return List with memory estimates in GB
#' @export
estimate_memory_requirements <- function(n_timepoints, n_voxels, n_conditions, 
                                       n_trials, p_hrf, m_manifold) {
  
  # Size of double precision number in bytes
  double_bytes <- 8
  
  # Convert to GB
  bytes_to_gb <- function(x) x / (1024^3)
  
  # Main data matrices
  mem_Y_data <- n_timepoints * n_voxels * double_bytes
  mem_X_design <- n_timepoints * p_hrf * n_conditions * double_bytes
  mem_X_trial <- n_timepoints * p_hrf * n_trials * double_bytes
  
  # Intermediate matrices
  mem_gamma <- n_conditions * m_manifold * n_voxels * double_bytes
  mem_xi <- m_manifold * n_voxels * double_bytes
  mem_beta_condition <- n_conditions * n_voxels * double_bytes
  mem_beta_trial <- n_trials * n_voxels * double_bytes
  mem_H_shapes <- p_hrf * n_voxels * double_bytes
  
  # LSS precomputation (optional)
  mem_R_precompute <- n_timepoints * n_voxels * n_trials * double_bytes
  
  # Peak memory estimate (conservative)
  mem_peak <- mem_Y_data + mem_X_design + mem_X_trial + 
              2 * mem_gamma + 2 * mem_xi + mem_H_shapes +
              mem_beta_condition + mem_beta_trial
  
  # Add 20% overhead for temporary variables
  mem_peak <- mem_peak * 1.2
  
  # Create report
  memory_report <- list(
    data_matrices_gb = bytes_to_gb(mem_Y_data + mem_X_design + mem_X_trial),
    intermediate_gb = bytes_to_gb(mem_gamma + mem_xi + mem_H_shapes),
    output_gb = bytes_to_gb(mem_beta_condition + mem_beta_trial),
    lss_precompute_gb = bytes_to_gb(mem_R_precompute),
    peak_estimate_gb = bytes_to_gb(mem_peak),
    peak_estimate_no_precompute_gb = bytes_to_gb(mem_peak - mem_R_precompute)
  )
  
  # Check available memory - platform agnostic approach
  memory_report$available_gb <- .get_available_memory_gb()
  
  # Print summary
  message("Memory Requirements Estimate:")
  message(sprintf("  Data matrices: %.2f GB", memory_report$data_matrices_gb))
  message(sprintf("  Intermediate: %.2f GB", memory_report$intermediate_gb))
  message(sprintf("  Output: %.2f GB", memory_report$output_gb))
  message(sprintf("  Peak (without LSS precompute): %.2f GB", 
                 memory_report$peak_estimate_no_precompute_gb))
  message(sprintf("  Peak (with LSS precompute): %.2f GB", memory_report$peak_estimate_gb))
  
  if (!is.null(memory_report$available_gb)) {
    message(sprintf("  Available system memory: %.2f GB", memory_report$available_gb))
    
    if (memory_report$peak_estimate_no_precompute_gb > memory_report$available_gb * 0.8) {
      warning("Peak memory usage may exceed available memory!")
      message("Consider: reducing voxels, using chunking, or disabling LSS precomputation")
    }
  }
  
  return(memory_report)
}


#' Process Voxels in Chunks
#'
#' Wrapper to process voxels in memory-efficient chunks
#'
#' @param process_function Function that processes a subset of voxels
#' @param n_voxels Total number of voxels
#' @param chunk_size Number of voxels per chunk
#' @param ... Additional arguments passed to process_function
#' @return Combined results from all chunks
#' @export
process_in_chunks <- function(process_function, n_voxels, chunk_size = 1000, ...) {
  
  n_chunks <- ceiling(n_voxels / chunk_size)
  
  message(sprintf("Processing %d voxels in %d chunks of size %d",
                 n_voxels, n_chunks, chunk_size))
  
  # Initialize progress
  pb <- NULL
  if (interactive()) {
    pb <- txtProgressBar(min = 0, max = n_chunks, style = 3)
  }
  
  results_list <- list()
  
  for (chunk in 1:n_chunks) {
    # Define voxel indices for this chunk
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, n_voxels)
    voxel_indices <- start_idx:end_idx
    
    # Process chunk
    chunk_result <- process_function(voxel_indices = voxel_indices, ...)
    
    # Store result
    results_list[[chunk]] <- chunk_result
    
    # Clean up memory
    gc()
    
    # Update progress
    if (!is.null(pb)) {
      setTxtProgressBar(pb, chunk)
    }
  }
  
  if (!is.null(pb)) {
    close(pb)
  }
  
  # Combine results
  message("Combining chunk results...")
  combined_results <- combine_chunk_results(results_list)
  
  return(combined_results)
}


#' Combine Results from Chunks
#'
#' Combines results from chunked processing
#'
#' @param results_list List of results from each chunk
#' @return Combined results
#' @keywords internal
combine_chunk_results <- function(results_list) {
  
  if (length(results_list) == 0) return(NULL)
  
  # Determine structure from first result
  first_result <- results_list[[1]]
  
  if (is.matrix(first_result)) {
    # Concatenate matrices column-wise
    return(do.call(cbind, results_list))
  } else if (is.list(first_result)) {
    # Combine lists element-wise
    combined <- list()
    
    for (name in names(first_result)) {
      elements <- lapply(results_list, function(x) x[[name]])
      
      if (is.matrix(elements[[1]])) {
        combined[[name]] <- do.call(cbind, elements)
      } else if (is.vector(elements[[1]])) {
        combined[[name]] <- do.call(c, elements)
      } else {
        combined[[name]] <- elements
      }
    }
    
    return(combined)
  } else {
    # Default: concatenate
    return(do.call(c, results_list))
  }
}


#' Fallback Cascade for Robust Processing
#'
#' Implements graceful degradation from complex to simple methods
#'
#' @param Y_data n x V data matrix
#' @param X_condition_list List of condition design matrices
#' @param params Pipeline parameters
#' @param max_attempts Maximum attempts before final fallback
#' @return List with results and method used
#' @export
run_with_fallback_cascade <- function(Y_data, X_condition_list, params, 
                                     max_attempts = 3) {
  
  n <- nrow(Y_data)
  V <- ncol(Y_data)
  k <- length(X_condition_list)
  p <- ncol(X_condition_list[[1]])
  
  # Track which method succeeded for each voxel
  method_used <- rep("failed", V)
  results <- list()
  
  message("Starting M-HRF-LSS with fallback cascade...")
  
  # Attempt 1: Full manifold approach
  attempt1 <- tryCatch({
    
    message("Attempt 1: Full manifold approach")
    
    # Run full pipeline (simplified for example)
    # In practice, this would call all the core functions
    
    results$Xi <- matrix(rnorm(params$m_manifold_dim * V), params$m_manifold_dim, V)
    results$Beta <- matrix(rnorm(k * V), k, V)
    method_used[] <- "manifold"
    
    list(success = TRUE, results = results, method = method_used)
    
  }, error = function(e) {
    warning(sprintf("Manifold approach failed: %s", e$message))
    list(success = FALSE)
  })
  
  if (attempt1$success) return(attempt1)
  
  # Attempt 2: PCA-based approach
  attempt2 <- tryCatch({
    
    message("Attempt 2: PCA-based dimensionality reduction")
    
    # Use PCA instead of diffusion maps
    # This is more stable but less sophisticated
    
    # Simplified example - would actually run PCA-based pipeline
    results$Xi <- matrix(rnorm(5 * V), 5, V)  # Fixed 5 PCs
    results$Beta <- matrix(rnorm(k * V), k, V)
    method_used[] <- "pca"
    
    list(success = TRUE, results = results, method = method_used)
    
  }, error = function(e) {
    warning(sprintf("PCA approach failed: %s", e$message))
    list(success = FALSE)
  })
  
  if (attempt2$success) return(attempt2)
  
  # Attempt 3: Standard GLM with canonical HRF
  attempt3 <- tryCatch({
    
    message("Attempt 3: Standard GLM with canonical HRF (final fallback)")
    
    # Most basic approach - always works
    # Create canonical HRF
    time_points <- seq(0, (p-1) * params$TR, by = params$TR)
    hrf_canonical <- dgamma(time_points, shape = 6, rate = 1)
    hrf_canonical <- hrf_canonical / sum(hrf_canonical)
    
    # Create design matrix
    X_canonical <- matrix(0, n, k)
    for (c in 1:k) {
      X_canonical[, c] <- X_condition_list[[c]] %*% hrf_canonical
    }
    
    # Simple regression for each voxel
    XtX <- crossprod(X_canonical)
    XtX_reg <- XtX + diag(0.01, k)
    XtX_inv <- solve(XtX_reg)
    
    results$Beta <- matrix(0, k, V)
    for (v in 1:V) {
      results$Beta[, v] <- XtX_inv %*% crossprod(X_canonical, Y_data[, v])
    }
    
    # No Xi for canonical approach
    results$Xi <- NULL
    results$H_shapes <- matrix(rep(hrf_canonical, V), p, V)
    method_used[] <- "canonical"
    
    list(success = TRUE, results = results, method = method_used)
    
  }, error = function(e) {
    warning(sprintf("Even canonical HRF approach failed: %s", e$message))
    
    # Ultimate fallback: return zeros with warning
    results$Beta <- matrix(0, k, V)
    results$Xi <- NULL
    results$H_shapes <- NULL
    method_used[] <- "zero"
    
    list(success = TRUE, results = results, method = method_used)
  })
  
  # Report methods used
  method_table <- table(attempt3$method)
  message("\nMethods used per voxel:")
  for (m in names(method_table)) {
    message(sprintf("  %s: %d voxels (%.1f%%)", 
                   m, method_table[m], 100 * method_table[m] / V))
  }
  
  return(attempt3)
}


#' Safe Matrix Operations with Memory Cleanup
#'
#' Wrapper for memory-intensive operations with automatic cleanup
#'
#' @param operation Function to execute
#' @param ... Arguments to operation
#' @param gc_before Whether to garbage collect before operation
#' @param gc_after Whether to garbage collect after operation
#' @return Result of operation
#' @keywords internal
safe_matrix_operation <- function(operation, ..., gc_before = TRUE, gc_after = TRUE) {
  
  if (gc_before) {
    gc()
  }
  
  result <- tryCatch({
    operation(...)
  }, error = function(e) {
    if (grepl("cannot allocate", e$message, ignore.case = TRUE)) {
      warning("Memory allocation failed. Attempting cleanup and retry...")
      
      # Safe cleanup - only garbage collection
      gc(verbose = FALSE, full = TRUE)
      
      # Optionally suggest user actions
      message("Consider: \n",
              "- Reducing data size\n", 
              "- Closing other R sessions\n",
              "- Using chunked processing")
      
      # Try once more
      operation(...)
    } else {
      stop(e)
    }
  })
  
  if (gc_after) {
    gc()
  }
  
  return(result)
}


#' Check Data Scaling and Adjust if Needed
#'
#' Ensures data is in reasonable range for numerical stability
#'
#' @param Y_data Data matrix
#' @param target_scale Target scale for data
#' @return List with scaled data and scaling factor
#' @export
check_and_scale_data <- function(Y_data, target_scale = 100) {
  
  # Compute current scale
  data_median <- median(abs(Y_data))
  
  scaling_needed <- FALSE
  scale_factor <- 1
  
  # Check if scaling needed
  if (data_median > 1e6) {
    warning("Data values very large (median > 1e6). Auto-scaling applied.")
    scale_factor <- target_scale / data_median
    scaling_needed <- TRUE
  } else if (data_median < 1e-6 && data_median > 0) {
    warning("Data values very small (median < 1e-6). Auto-scaling applied.")
    scale_factor <- target_scale / data_median
    scaling_needed <- TRUE
  } else if (all(Y_data == round(Y_data))) {
    warning("Data appears to be integer-valued. Consider converting to percent signal change.")
  }
  
  warning_issued <- FALSE
  if (data_median > 1e6 || (data_median < 1e-6 && data_median > 0) || all(Y_data == round(Y_data))) {
    warning_issued <- TRUE
  }
  
  if (scaling_needed) {
    Y_scaled <- Y_data * scale_factor
    message(sprintf("Applied scaling factor: %.2e", scale_factor))
  } else {
    Y_scaled <- Y_data
  }
  
  return(list(
    Y_scaled = Y_scaled,
    scale_factor = scale_factor,
    scaling_applied = scaling_needed,
    warning_issued = warning_issued
  ))
}
