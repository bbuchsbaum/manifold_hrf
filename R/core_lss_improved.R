# Improved Core LSS Implementation with Memory Management

#' Run LSS Voxel Loop with Intelligent Memory Management (Core)
#'
#' Enhanced version of run_lss_voxel_loop_core that actually uses
#' the ram_heuristic_GB_for_Rt parameter to decide computation strategy
#'
#' @inheritParams run_lss_voxel_loop_core
#' @return A T x V matrix of trial-wise beta estimates
#' @export
run_lss_voxel_loop_core_improved <- function(Y_proj_matrix,
                                           X_trial_onset_list_of_matrices,
                                           H_shapes_allvox_matrix,
                                           A_lss_fixed_matrix,
                                           P_lss_matrix,
                                           p_lss_vector,
                                           n_jobs = 1,
                                           ram_heuristic_GB_for_Rt = 1.0) {
  
  if (!is.matrix(Y_proj_matrix)) {
    stop("Y_proj_matrix must be a matrix")
  }
  
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  # Calculate memory requirement for full precomputation
  bytes_per_double <- 8
  precompute_memory_gb <- (n * T_trials * V * bytes_per_double) / (1024^3)
  
  # Log memory analysis
  if (getOption("manifoldhrf.verbose", TRUE)) {
    message(sprintf(
      "\nLSS Memory Analysis:",
      "\n  Data dimensions: %d timepoints × %d trials × %d voxels",
      n, T_trials, V
    ))
    message(sprintf(
      "  Full precomputation would require: %.2f GB",
      precompute_memory_gb
    ))
    message(sprintf(
      "  Memory limit (ram_heuristic_GB_for_Rt): %.2f GB",
      ram_heuristic_GB_for_Rt
    ))
  }
  
  # STRATEGY 1: Full precomputation (current implementation)
  if (precompute_memory_gb <= ram_heuristic_GB_for_Rt) {
    if (getOption("manifoldhrf.verbose", TRUE)) {
      message("  → Strategy: FULL PRECOMPUTATION (fastest)")
    }
    
    # Pre-compute all convolved trial regressors (n x T x V tensor)
    all_C <- array(0, dim = c(n, T_trials, V))
    for (t in seq_len(T_trials)) {
      all_C[, t, ] <- X_trial_onset_list_of_matrices[[t]] %*% H_shapes_allvox_matrix
    }
    
    # Process voxels in parallel
    voxel_fun <- function(v) {
      tryCatch({
        C_v <- all_C[, , v]  # n x T matrix
        result <- fmrilss::lss(
          Y = Y_proj_matrix[, v, drop = FALSE],
          X = C_v,
          Z = NULL,  # Data is already projected
          method = "r_optimized"
        )
        as.vector(result)
      }, error = function(e) {
        warning(sprintf("LSS failed for voxel %d: %s", v, e$message))
        rep(NA_real_, T_trials)
      })
    }
    
    res_list <- .parallel_lapply(seq_len(V), voxel_fun, n_jobs)
    return(do.call(cbind, res_list))
  }
  
  # STRATEGY 2: Chunked processing
  # Calculate optimal chunk size
  memory_per_voxel_gb <- precompute_memory_gb / V
  optimal_chunk_size <- floor(ram_heuristic_GB_for_Rt / memory_per_voxel_gb)
  
  if (optimal_chunk_size >= 100) {  # Chunking is viable
    optimal_chunk_size <- min(optimal_chunk_size, 1000)  # Cap at 1000 for efficiency
    n_chunks <- ceiling(V / optimal_chunk_size)
    
    if (getOption("manifoldhrf.verbose", TRUE)) {
      message(sprintf(
        "  → Strategy: CHUNKED PROCESSING (%d chunks of ~%d voxels)",
        n_chunks, optimal_chunk_size
      ))
    }
    
    # Initialize result matrix
    Beta_all <- matrix(NA_real_, T_trials, V)
    
    # Process chunks
    if (interactive() && getOption("manifoldhrf.progress", TRUE)) {
      pb <- txtProgressBar(min = 0, max = n_chunks, style = 3)
    } else {
      pb <- NULL
    }
    
    for (chunk_idx in seq_len(n_chunks)) {
      # Define voxel range
      start_v <- (chunk_idx - 1) * optimal_chunk_size + 1
      end_v <- min(chunk_idx * optimal_chunk_size, V)
      voxel_range <- start_v:end_v
      chunk_size <- length(voxel_range)
      
      # Precompute for this chunk only
      chunk_C <- array(0, dim = c(n, T_trials, chunk_size))
      H_chunk <- H_shapes_allvox_matrix[, voxel_range, drop = FALSE]
      
      for (t in seq_len(T_trials)) {
        chunk_C[, t, ] <- X_trial_onset_list_of_matrices[[t]] %*% H_chunk
      }
      
      # Process voxels in chunk
      chunk_fun <- function(v_local) {
        tryCatch({
          result <- fmrilss::lss(
            Y = Y_proj_matrix[, voxel_range[v_local], drop = FALSE],
            X = chunk_C[, , v_local],
            Z = NULL,
            method = "r_optimized"
          )
          as.vector(result)
        }, error = function(e) {
          warning(sprintf("LSS failed for voxel %d: %s", 
                         voxel_range[v_local], e$message))
          rep(NA_real_, T_trials)
        })
      }
      
      chunk_results <- .parallel_lapply(seq_len(chunk_size), chunk_fun, n_jobs)
      Beta_all[, voxel_range] <- do.call(cbind, chunk_results)
      
      # Memory cleanup
      rm(chunk_C, H_chunk)
      if (chunk_idx %% 5 == 0) gc()  # Periodic garbage collection
      
      if (!is.null(pb)) setTxtProgressBar(pb, chunk_idx)
    }
    
    if (!is.null(pb)) close(pb)
    return(Beta_all)
  }
  
  # STRATEGY 3: Pure streaming (no precomputation)
  if (getOption("manifoldhrf.verbose", TRUE)) {
    message("  → Strategy: STREAMING (most memory-efficient)")
  }
  
  voxel_fun <- function(v) {
    tryCatch({
      # Compute convolved regressors on-the-fly
      C_v <- matrix(0, n, T_trials)
      h_v <- H_shapes_allvox_matrix[, v]
      
      for (t in seq_len(T_trials)) {
        C_v[, t] <- X_trial_onset_list_of_matrices[[t]] %*% h_v
      }
      
      result <- fmrilss::lss(
        Y = Y_proj_matrix[, v, drop = FALSE],
        X = C_v,
        Z = NULL,
        method = "r_optimized"
      )
      as.vector(result)
    }, error = function(e) {
      warning(sprintf("LSS failed for voxel %d: %s", v, e$message))
      rep(NA_real_, T_trials)
    })
  }
  
  res_list <- .parallel_lapply(seq_len(V), voxel_fun, n_jobs)
  do.call(cbind, res_list)
}

#' Calculate Optimal ram_heuristic_GB_for_Rt
#'
#' Provides intelligent default based on system resources and data size
#'
#' @param n_timepoints Number of timepoints
#' @param n_voxels Number of voxels  
#' @param n_trials Number of trials
#' @param target_strategy One of "fast", "balanced", "memory" 
#' @param safety_margin Fraction of available RAM to use (0-1)
#' @return Recommended ram_heuristic_GB_for_Rt value
#' @export
calculate_optimal_ram_heuristic <- function(n_timepoints, 
                                          n_voxels,
                                          n_trials,
                                          target_strategy = "balanced",
                                          safety_margin = 0.5) {
  
  # Get system memory
  sys_info <- Sys.info()
  
  if (sys_info["sysname"] == "Darwin") {  # macOS
    total_ram_bytes <- as.numeric(system("sysctl -n hw.memsize", intern = TRUE))
  } else if (sys_info["sysname"] == "Linux") {
    meminfo <- readLines("/proc/meminfo")
    mem_line <- grep("^MemTotal:", meminfo, value = TRUE)
    total_ram_kb <- as.numeric(gsub("[^0-9]", "", mem_line))
    total_ram_bytes <- total_ram_kb * 1024
  } else {  # Windows
    total_ram_mb <- memory.limit()
    total_ram_bytes <- total_ram_mb * 1024 * 1024
  }
  
  total_ram_gb <- total_ram_bytes / (1024^3)
  available_ram_gb <- total_ram_gb * safety_margin
  
  # Calculate memory for different strategies
  bytes_per_double <- 8
  full_precompute_gb <- (n_timepoints * n_trials * n_voxels * bytes_per_double) / (1024^3)
  
  # Strategy-specific recommendations
  if (target_strategy == "fast") {
    # Prefer full precomputation if it fits
    recommended <- ifelse(full_precompute_gb < available_ram_gb,
                         full_precompute_gb * 1.2,  # 20% overhead
                         available_ram_gb * 0.8)     # 80% of available
  } else if (target_strategy == "memory") {
    # Conservative: force chunking
    recommended <- min(1.0, available_ram_gb * 0.2)
  } else {  # balanced
    # Aim for reasonable chunk sizes
    target_chunk_size <- 500  # voxels
    chunk_memory_gb <- (n_timepoints * n_trials * target_chunk_size * bytes_per_double) / (1024^3)
    recommended <- min(chunk_memory_gb * 1.5, available_ram_gb * 0.5)
  }
  
  message(sprintf(
    "\nRAM Heuristic Calculation:",
    "\n  System RAM: %.1f GB total, %.1f GB available (%.0f%% margin)",
    total_ram_gb, available_ram_gb, safety_margin * 100
  ))
  message(sprintf(
    "  Full precomputation needs: %.2f GB",
    full_precompute_gb
  ))
  message(sprintf(
    "  Strategy '%s' recommends: %.2f GB",
    target_strategy, recommended
  ))
  
  return(recommended)
}

#' Memory-Efficient Data Structure for Trial Regressors
#'
#' Alternative to 3D array using sparse representations
#'
#' @keywords internal
create_sparse_trial_regressors <- function(X_trial_list, H_matrix) {
  n <- nrow(X_trial_list[[1]])
  p <- ncol(X_trial_list[[1]])
  T_trials <- length(X_trial_list)
  V <- ncol(H_matrix)
  
  # Use Matrix package for sparse storage
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package required for sparse storage")
  }
  
  # Create list of sparse matrices instead of 3D array
  sparse_C_list <- vector("list", T_trials)
  
  for (t in seq_len(T_trials)) {
    # X_trial is often sparse (single trial onset)
    X_sparse <- Matrix::Matrix(X_trial_list[[t]], sparse = TRUE)
    
    # Compute C for all voxels for this trial
    C_t <- X_sparse %*% H_matrix  # n x V
    
    # Store as sparse if beneficial
    sparsity <- sum(C_t == 0) / length(C_t)
    if (sparsity > 0.5) {
      sparse_C_list[[t]] <- Matrix::Matrix(C_t, sparse = TRUE)
    } else {
      sparse_C_list[[t]] <- C_t
    }
  }
  
  # Return accessor function instead of full array
  get_C_for_voxel <- function(v) {
    C_v <- matrix(0, n, T_trials)
    for (t in seq_len(T_trials)) {
      C_v[, t] <- sparse_C_list[[t]][, v]
    }
    return(C_v)
  }
  
  structure(
    list(
      get_voxel = get_C_for_voxel,
      sparse_list = sparse_C_list
    ),
    class = "sparse_trial_regressors"
  )
}