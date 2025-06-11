# LSS Memory Optimization Examples
# Concrete implementations for intelligent memory management

#' Run LSS with Intelligent Memory Management
#'
#' Modified version of run_lss_voxel_loop_core that uses ram_heuristic_GB_for_Rt
#' to decide between precomputation and streaming computation
#'
#' @inheritParams run_lss_voxel_loop_core
#' @export
run_lss_voxel_loop_memory_aware <- function(Y_proj_matrix,
                                           X_trial_onset_list_of_matrices,
                                           H_shapes_allvox_matrix,
                                           A_lss_fixed_matrix,
                                           P_lss_matrix,
                                           p_lss_vector,
                                           n_jobs = 1,
                                           ram_heuristic_GB_for_Rt = 1.0) {
  
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  # Calculate memory requirement for precomputation
  bytes_per_double <- 8
  precompute_memory_gb <- (n * T_trials * V * bytes_per_double) / (1024^3)
  
  message(sprintf(
    "LSS memory analysis: %d timepoints × %d trials × %d voxels",
    n, T_trials, V
  ))
  message(sprintf(
    "  Precomputation would require: %.2f GB",
    precompute_memory_gb
  ))
  message(sprintf(
    "  RAM heuristic limit: %.2f GB",
    ram_heuristic_GB_for_Rt
  ))
  
  # Decision: precompute or stream?
  use_precompute <- precompute_memory_gb <= ram_heuristic_GB_for_Rt
  
  if (use_precompute) {
    message("  → Using precomputation strategy (fits in memory)")
    return(run_lss_with_precomputation(
      Y_proj_matrix, X_trial_onset_list_of_matrices, H_shapes_allvox_matrix,
      A_lss_fixed_matrix, P_lss_matrix, p_lss_vector, n_jobs
    ))
  } else {
    message("  → Using streaming strategy (memory-efficient)")
    return(run_lss_with_streaming(
      Y_proj_matrix, X_trial_onset_list_of_matrices, H_shapes_allvox_matrix,
      A_lss_fixed_matrix, P_lss_matrix, p_lss_vector, n_jobs
    ))
  }
}

#' Run LSS with Precomputation (Current Implementation)
#' @keywords internal
run_lss_with_precomputation <- function(Y_proj_matrix,
                                       X_trial_onset_list_of_matrices,
                                       H_shapes_allvox_matrix,
                                       A_lss_fixed_matrix,
                                       P_lss_matrix,
                                       p_lss_vector,
                                       n_jobs) {
  
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  # Precompute all convolved trial regressors
  all_C <- array(0, dim = c(n, T_trials, V))
  for (t in seq_len(T_trials)) {
    all_C[, t, ] <- X_trial_onset_list_of_matrices[[t]] %*% H_shapes_allvox_matrix
  }
  
  # Process voxels
  voxel_fun <- function(v) {
    C_v <- all_C[, , v]
    fmrilss::lss(
      Y = Y_proj_matrix[, v, drop = FALSE],
      X = C_v,
      Z = NULL,
      method = "r_optimized"
    )
  }
  
  res_list <- .parallel_lapply(seq_len(V), voxel_fun, n_jobs)
  do.call(cbind, res_list)
}

#' Run LSS with Streaming (Memory-Efficient)
#' @keywords internal
run_lss_with_streaming <- function(Y_proj_matrix,
                                 X_trial_onset_list_of_matrices,
                                 H_shapes_allvox_matrix,
                                 A_lss_fixed_matrix,
                                 P_lss_matrix,
                                 p_lss_vector,
                                 n_jobs) {
  
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  # Process voxels without precomputation
  voxel_fun <- function(v) {
    # Compute convolved regressors on-the-fly for this voxel only
    C_v <- matrix(0, n, T_trials)
    h_v <- H_shapes_allvox_matrix[, v]
    
    for (t in seq_len(T_trials)) {
      C_v[, t] <- X_trial_onset_list_of_matrices[[t]] %*% h_v
    }
    
    fmrilss::lss(
      Y = Y_proj_matrix[, v, drop = FALSE],
      X = C_v,
      Z = NULL,
      method = "r_optimized"
    )
  }
  
  res_list <- .parallel_lapply(seq_len(V), voxel_fun, n_jobs)
  do.call(cbind, res_list)
}

#' Intelligent RAM Heuristic Calculator
#'
#' Calculates appropriate ram_heuristic_GB_for_Rt based on system resources
#' and dataset characteristics
#'
#' @param safety_factor Fraction of available RAM to use (default 0.5)
#' @return Recommended ram_heuristic_GB_for_Rt value
#' @export
calculate_ram_heuristic <- function(safety_factor = 0.5) {
  # Get available memory
  if (Sys.info()["sysname"] == "Windows") {
    # Windows: use memory.limit()
    available_mb <- memory.limit()
    available_gb <- available_mb / 1024
  } else if (Sys.info()["sysname"] == "Darwin") {
    # macOS: use system command
    total_gb <- as.numeric(system("sysctl -n hw.memsize", intern = TRUE)) / (1024^3)
    # Estimate available as 50% of total (conservative)
    available_gb <- total_gb * 0.5
  } else {
    # Linux: parse /proc/meminfo
    meminfo <- readLines("/proc/meminfo")
    available_line <- grep("^MemAvailable:", meminfo, value = TRUE)
    if (length(available_line) > 0) {
      available_kb <- as.numeric(gsub("[^0-9]", "", available_line))
      available_gb <- available_kb / (1024^2)
    } else {
      # Fallback: use total memory
      total_line <- grep("^MemTotal:", meminfo, value = TRUE)
      total_kb <- as.numeric(gsub("[^0-9]", "", total_line))
      available_gb <- total_kb / (1024^2) * 0.5
    }
  }
  
  # Apply safety factor
  recommended_gb <- available_gb * safety_factor
  
  message(sprintf(
    "System memory analysis:\n  Available: %.1f GB\n  Recommended limit: %.1f GB (%.0f%% safety factor)",
    available_gb, recommended_gb, safety_factor * 100
  ))
  
  return(recommended_gb)
}

#' Memory-Efficient LSS using Chunked Voxel Processing
#'
#' Process voxels in chunks to limit peak memory usage
#'
#' @inheritParams run_lss_voxel_loop_core
#' @param voxel_chunk_size Number of voxels to process at once
#' @export
run_lss_voxel_loop_chunked <- function(Y_proj_matrix,
                                      X_trial_onset_list_of_matrices,
                                      H_shapes_allvox_matrix,
                                      A_lss_fixed_matrix,
                                      P_lss_matrix,
                                      p_lss_vector,
                                      n_jobs = 1,
                                      voxel_chunk_size = 1000) {
  
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  # Calculate optimal chunk size based on memory
  bytes_per_element <- 8
  memory_per_voxel_mb <- (n * T_trials * bytes_per_element) / (1024^2)
  available_mb <- 1024  # 1 GB target per chunk
  optimal_chunk_size <- floor(available_mb / memory_per_voxel_mb)
  
  # Use the smaller of user-specified and optimal
  chunk_size <- min(voxel_chunk_size, optimal_chunk_size, V)
  n_chunks <- ceiling(V / chunk_size)
  
  message(sprintf(
    "Chunked LSS: %d voxels in %d chunks of size %d",
    V, n_chunks, chunk_size
  ))
  
  # Initialize result matrix
  Beta_all <- matrix(NA_real_, T_trials, V)
  
  # Process chunks
  pb <- txtProgressBar(min = 0, max = n_chunks, style = 3)
  
  for (chunk_idx in seq_len(n_chunks)) {
    # Define voxel range for this chunk
    start_v <- (chunk_idx - 1) * chunk_size + 1
    end_v <- min(chunk_idx * chunk_size, V)
    voxel_range <- start_v:end_v
    chunk_V <- length(voxel_range)
    
    # Extract data for this chunk
    Y_chunk <- Y_proj_matrix[, voxel_range, drop = FALSE]
    H_chunk <- H_shapes_allvox_matrix[, voxel_range, drop = FALSE]
    
    # Precompute convolved regressors for this chunk only
    chunk_C <- array(0, dim = c(n, T_trials, chunk_V))
    for (t in seq_len(T_trials)) {
      chunk_C[, t, ] <- X_trial_onset_list_of_matrices[[t]] %*% H_chunk
    }
    
    # Process voxels in chunk
    chunk_fun <- function(v_local) {
      fmrilss::lss(
        Y = Y_chunk[, v_local, drop = FALSE],
        X = chunk_C[, , v_local],
        Z = NULL,
        method = "r_optimized"
      )
    }
    
    chunk_results <- .parallel_lapply(seq_len(chunk_V), chunk_fun, n_jobs)
    Beta_all[, voxel_range] <- do.call(cbind, chunk_results)
    
    # Clean up memory
    rm(chunk_C, Y_chunk, H_chunk)
    gc()
    
    setTxtProgressBar(pb, chunk_idx)
  }
  
  close(pb)
  return(Beta_all)
}

#' Future-Based Memory-Aware LSS
#'
#' Uses future package's memory limits for automatic resource management
#'
#' @inheritParams run_lss_voxel_loop_core
#' @param memory_limit_gb Memory limit per worker in GB
#' @export
run_lss_with_future_memory_limits <- function(Y_proj_matrix,
                                             X_trial_onset_list_of_matrices,
                                             H_shapes_allvox_matrix,
                                             A_lss_fixed_matrix,
                                             P_lss_matrix,
                                             p_lss_vector,
                                             n_jobs = 1,
                                             memory_limit_gb = 2) {
  
  # Check if future is available
  if (!requireNamespace("future", quietly = TRUE)) {
    warning("future package not available, falling back to standard parallel processing")
    return(run_lss_voxel_loop_core(
      Y_proj_matrix, X_trial_onset_list_of_matrices, H_shapes_allvox_matrix,
      A_lss_fixed_matrix, P_lss_matrix, p_lss_vector, n_jobs
    ))
  }
  
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  # Set up future with memory constraints
  old_plan <- future::plan()
  on.exit(future::plan(old_plan))
  
  # Configure workers with memory limits
  future::plan(
    future::multisession,
    workers = n_jobs,
    globals = list(
      maxSize = memory_limit_gb * 1024^3  # Convert GB to bytes
    )
  )
  
  # Create futures for voxel processing
  voxel_futures <- list()
  
  for (v in seq_len(V)) {
    # Extract voxel-specific data
    y_v <- Y_proj_matrix[, v]
    h_v <- H_shapes_allvox_matrix[, v]
    
    # Create future for this voxel
    voxel_futures[[v]] <- future::future({
      # Compute convolved regressors
      C_v <- matrix(0, n, T_trials)
      for (t in seq_len(T_trials)) {
        C_v[, t] <- X_trial_onset_list_of_matrices[[t]] %*% h_v
      }
      
      # Run LSS
      fmrilss::lss(
        Y = matrix(y_v, ncol = 1),
        X = C_v,
        Z = NULL,
        method = "r_optimized"
      )
    })
  }
  
  # Collect results
  Beta_all <- matrix(NA_real_, T_trials, V)
  
  pb <- txtProgressBar(min = 0, max = V, style = 3)
  for (v in seq_len(V)) {
    Beta_all[, v] <- future::value(voxel_futures[[v]])
    setTxtProgressBar(pb, v)
  }
  close(pb)
  
  return(Beta_all)
}

#' Hybrid Memory Strategy for LSS
#'
#' Combines precomputation for some voxels with streaming for others
#' based on available memory
#'
#' @inheritParams run_lss_voxel_loop_core
#' @export
run_lss_hybrid_memory <- function(Y_proj_matrix,
                                 X_trial_onset_list_of_matrices,
                                 H_shapes_allvox_matrix,
                                 A_lss_fixed_matrix,
                                 P_lss_matrix,
                                 p_lss_vector,
                                 n_jobs = 1,
                                 ram_heuristic_GB_for_Rt = 1.0) {
  
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  # Calculate memory per voxel
  bytes_per_double <- 8
  memory_per_voxel_gb <- (n * T_trials * bytes_per_double) / (1024^3)
  
  # How many voxels can we precompute?
  max_precompute_voxels <- floor(ram_heuristic_GB_for_Rt / memory_per_voxel_gb)
  max_precompute_voxels <- min(max_precompute_voxels, V)
  
  message(sprintf(
    "Hybrid strategy: precompute %d/%d voxels (%.1f%%), stream the rest",
    max_precompute_voxels, V, 100 * max_precompute_voxels / V
  ))
  
  # Initialize result
  Beta_all <- matrix(NA_real_, T_trials, V)
  
  if (max_precompute_voxels > 0) {
    # Phase 1: Precomputed voxels
    message("Phase 1: Processing precomputed voxels...")
    
    # Precompute convolved regressors for first set of voxels
    precompute_C <- array(0, dim = c(n, T_trials, max_precompute_voxels))
    H_subset <- H_shapes_allvox_matrix[, 1:max_precompute_voxels, drop = FALSE]
    
    for (t in seq_len(T_trials)) {
      precompute_C[, t, ] <- X_trial_onset_list_of_matrices[[t]] %*% H_subset
    }
    
    # Process precomputed voxels
    precompute_fun <- function(v) {
      fmrilss::lss(
        Y = Y_proj_matrix[, v, drop = FALSE],
        X = precompute_C[, , v],
        Z = NULL,
        method = "r_optimized"
      )
    }
    
    res_list <- .parallel_lapply(1:max_precompute_voxels, precompute_fun, n_jobs)
    Beta_all[, 1:max_precompute_voxels] <- do.call(cbind, res_list)
    
    # Clean up
    rm(precompute_C, H_subset)
    gc()
  }
  
  if (max_precompute_voxels < V) {
    # Phase 2: Streaming voxels
    remaining_voxels <- (max_precompute_voxels + 1):V
    message(sprintf("Phase 2: Processing %d streaming voxels...", length(remaining_voxels)))
    
    stream_fun <- function(v) {
      # Compute on-the-fly
      C_v <- matrix(0, n, T_trials)
      h_v <- H_shapes_allvox_matrix[, v]
      
      for (t in seq_len(T_trials)) {
        C_v[, t] <- X_trial_onset_list_of_matrices[[t]] %*% h_v
      }
      
      fmrilss::lss(
        Y = Y_proj_matrix[, v, drop = FALSE],
        X = C_v,
        Z = NULL,
        method = "r_optimized"
      )
    }
    
    res_list <- .parallel_lapply(remaining_voxels, stream_fun, n_jobs)
    Beta_all[, remaining_voxels] <- do.call(cbind, res_list)
  }
  
  return(Beta_all)
}