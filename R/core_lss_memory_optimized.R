# Memory-Optimized LSS Implementation
# Based on Gemini's recommendations for intelligent memory management

#' Calculate Optimal RAM Heuristic
#'
#' Determines the appropriate RAM limit based on dataset size and available memory
#'
#' @param n_timepoints Number of timepoints
#' @param n_voxels Number of voxels
#' @param n_trials Number of trials
#' @param strategy Character: "fast", "balanced", or "memory"
#' @return Recommended RAM limit in GB
#' @export
calculate_optimal_ram_heuristic <- function(n_timepoints, n_voxels, n_trials,
                                          strategy = c("balanced", "fast", "memory")) {
  strategy <- match.arg(strategy)
  
  # Get available memory
  available_gb <- as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", 
                                   intern = TRUE)) / 1024 / 1024
  if (is.na(available_gb)) {
    # Fallback for non-Linux systems
    available_gb <- 8  # Conservative default
  }
  
  # Calculate memory needed for full precomputation
  mem_needed_gb <- (n_timepoints * n_trials * n_voxels * 8) / 1e9
  
  # Strategy-based recommendations
  ram_heuristic <- switch(strategy,
    fast = available_gb * 0.6,      # Use 60% of RAM for speed
    balanced = available_gb * 0.4,   # Use 40% for balance
    memory = available_gb * 0.25     # Use 25% for memory safety
  )
  
  # Provide informative message
  if (mem_needed_gb > ram_heuristic) {
    message(sprintf(
      "Dataset requires %.1f GB for full precomputation, but only %.1f GB allocated.\n",
      mem_needed_gb, ram_heuristic,
      "Will use chunked/streaming processing."
    ))
  }
  
  return(ram_heuristic)
}

#' Memory-Optimized LSS Voxel Loop
#'
#' Improved version that actually uses ram_heuristic_GB_for_Rt parameter
#'
#' @inheritParams run_lss_voxel_loop_core
#' @export
run_lss_voxel_loop_memory_optimized <- function(Y_proj_matrix,
                                               X_trial_onset_list_of_matrices,
                                               H_shapes_allvox_matrix,
                                               A_lss_fixed_matrix,
                                               P_lss_matrix,
                                               p_lss_vector,
                                               n_jobs = 1,
                                               ram_heuristic_GB_for_Rt = NULL) {
  
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  # Auto-calculate RAM heuristic if not provided
  if (is.null(ram_heuristic_GB_for_Rt)) {
    ram_heuristic_GB_for_Rt <- calculate_optimal_ram_heuristic(
      n, V, T_trials, strategy = "balanced"
    )
  }
  
  # Calculate memory requirement for full precomputation
  mem_required_gb <- (n * T_trials * V * 8) / 1e9
  
  # Choose strategy based on available memory
  if (mem_required_gb <= ram_heuristic_GB_for_Rt) {
    # Strategy 1: Full precomputation (current approach)
    message("Using full precomputation strategy")
    return(run_lss_full_precompute(Y_proj_matrix, X_trial_onset_list_of_matrices,
                                  H_shapes_allvox_matrix, n_jobs))
    
  } else if (mem_required_gb <= ram_heuristic_GB_for_Rt * 10) {
    # Strategy 2: Chunked processing
    chunk_size <- floor(V * ram_heuristic_GB_for_Rt / mem_required_gb)
    chunk_size <- max(100, min(chunk_size, 1000))  # Reasonable bounds
    message(sprintf("Using chunked processing with chunk size %d", chunk_size))
    return(run_lss_chunked(Y_proj_matrix, X_trial_onset_list_of_matrices,
                          H_shapes_allvox_matrix, chunk_size, n_jobs))
    
  } else {
    # Strategy 3: Streaming (no precomputation)
    message("Using streaming strategy (minimal memory usage)")
    return(run_lss_streaming(Y_proj_matrix, X_trial_onset_list_of_matrices,
                           H_shapes_allvox_matrix, n_jobs))
  }
}

#' Full Precomputation Strategy (Original)
#' @keywords internal
run_lss_full_precompute <- function(Y_proj_matrix, X_trial_onset_list_of_matrices,
                                   H_shapes_allvox_matrix, n_jobs) {
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  # Original approach - precompute everything
  all_C <- array(0, dim = c(n, T_trials, V))
  for (t in seq_len(T_trials)) {
    all_C[, t, ] <- X_trial_onset_list_of_matrices[[t]] %*% H_shapes_allvox_matrix
  }
  
  # Process voxels
  voxel_fun <- function(v) {
    C_v <- all_C[, , v]
    result <- fmrilss::lss(
      Y = Y_proj_matrix[, v, drop = FALSE],
      X = C_v,
      Z = NULL,
      method = "r_optimized"
    )
    as.vector(result)
  }
  
  res_list <- .parallel_lapply(seq_len(V), voxel_fun, n_jobs)
  do.call(cbind, res_list)
}

#' Chunked Processing Strategy
#' @keywords internal
run_lss_chunked <- function(Y_proj_matrix, X_trial_onset_list_of_matrices,
                           H_shapes_allvox_matrix, chunk_size, n_jobs) {
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  T_trials <- length(X_trial_onset_list_of_matrices)
  n_chunks <- ceiling(V / chunk_size)
  
  # Process in chunks
  Beta_matrix <- matrix(NA, T_trials, V)
  
  pb <- progress::progress_bar$new(
    format = "Processing chunks [:bar] :percent eta: :eta",
    total = n_chunks
  )
  
  for (chunk in seq_len(n_chunks)) {
    start_v <- (chunk - 1) * chunk_size + 1
    end_v <- min(chunk * chunk_size, V)
    chunk_voxels <- start_v:end_v
    chunk_size_actual <- length(chunk_voxels)
    
    # Precompute only for this chunk
    chunk_C <- array(0, dim = c(n, T_trials, chunk_size_actual))
    for (t in seq_len(T_trials)) {
      chunk_C[, t, ] <- X_trial_onset_list_of_matrices[[t]] %*% 
                        H_shapes_allvox_matrix[, chunk_voxels]
    }
    
    # Process chunk voxels
    chunk_fun <- function(v_idx) {
      C_v <- chunk_C[, , v_idx]
      v <- chunk_voxels[v_idx]
      result <- fmrilss::lss(
        Y = Y_proj_matrix[, v, drop = FALSE],
        X = C_v,
        Z = NULL,
        method = "r_optimized"
      )
      as.vector(result)
    }
    
    chunk_results <- .parallel_lapply(seq_len(chunk_size_actual), chunk_fun, n_jobs)
    Beta_matrix[, chunk_voxels] <- do.call(cbind, chunk_results)
    
    pb$tick()
  }
  
  return(Beta_matrix)
}

#' Streaming Strategy (Minimal Memory)
#' @keywords internal
run_lss_streaming <- function(Y_proj_matrix, X_trial_onset_list_of_matrices,
                             H_shapes_allvox_matrix, n_jobs) {
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  # Process each voxel without precomputation
  voxel_fun <- function(v) {
    # Compute C_v on the fly
    C_v <- matrix(0, n, T_trials)
    for (t in seq_len(T_trials)) {
      C_v[, t] <- X_trial_onset_list_of_matrices[[t]] %*% H_shapes_allvox_matrix[, v]
    }
    
    result <- fmrilss::lss(
      Y = Y_proj_matrix[, v, drop = FALSE],
      X = C_v,
      Z = NULL,
      method = "r_optimized"
    )
    as.vector(result)
  }
  
  # Add progress bar for long computations
  if (V > 1000) {
    voxel_fun <- progressr::with_progress({
      p <- progressr::progressor(V)
      function(v) {
        p()
        voxel_fun(v)
      }
    })
  }
  
  res_list <- .parallel_lapply(seq_len(V), voxel_fun, n_jobs)
  do.call(cbind, res_list)
}

#' Run LSS with Future Memory Limits
#'
#' Uses future package's resource management
#'
#' @export
run_lss_with_future_memory_limits <- function(Y_proj_matrix, 
                                             X_trial_onset_list_of_matrices,
                                             H_shapes_allvox_matrix,
                                             max_memory_gb = NULL) {
  if (is.null(max_memory_gb)) {
    max_memory_gb <- calculate_optimal_ram_heuristic(
      nrow(Y_proj_matrix),
      ncol(Y_proj_matrix),
      length(X_trial_onset_list_of_matrices),
      strategy = "balanced"
    )
  }
  
  # Configure future with memory limits
  oplan <- future::plan()
  on.exit(future::plan(oplan))
  
  future::plan(future::multisession, workers = 4, gc = TRUE,
               globals = structure(TRUE, 
                 maxSize = max_memory_gb * 1024^3))
  
  # Use future.apply for memory-aware parallelization
  V <- ncol(Y_proj_matrix)
  
  # Let future handle chunking based on memory constraints
  Beta_matrix <- future.apply::future_sapply(seq_len(V), function(v) {
    n <- nrow(Y_proj_matrix)
    T_trials <- length(X_trial_onset_list_of_matrices)
    
    # Compute on demand
    C_v <- matrix(0, n, T_trials)
    for (t in seq_len(T_trials)) {
      C_v[, t] <- X_trial_onset_list_of_matrices[[t]] %*% H_shapes_allvox_matrix[, v]
    }
    
    result <- fmrilss::lss(
      Y = Y_proj_matrix[, v, drop = FALSE],
      X = C_v,
      Z = NULL,
      method = "r_optimized"
    )
    as.vector(result)
  }, future.seed = TRUE)
  
  return(Beta_matrix)
}