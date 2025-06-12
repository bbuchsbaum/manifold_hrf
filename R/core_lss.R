# Core LSS (Least Squares Separate) Functions - Consolidated Version
# This file consolidates all LSS functionality from the previous 4 files
# into a single, well-organized implementation with multiple memory strategies

#' Run LSS voxel loop with memory optimization strategies
#'
#' Main dispatcher function that performs trial-wise LSS estimation using
#' one of three memory strategies based on dataset size and available resources.
#'
#' @param Y_proj_matrix n x V matrix of projected BOLD data
#' @param X_trial_onset_list_of_matrices List of T sparse matrices (n x p each)
#' @param B_reconstructor_matrix p x m manifold basis matrix
#' @param Xi_smoothed_allvox_matrix m x V smoothed manifold coordinates
#' @param A_lss_fixed_matrix Optional n x q fixed regressors matrix
#' @param memory_strategy Character string specifying computation method:
#'   \itemize{
#'     \item{"auto"}{ (Default) Automatically selects based on data size}
#'     \item{"full"}{ Precomputes all trial regressors. Fastest but uses most memory}
#'     \item{"chunked"}{ Processes trials in batches. Good balance of speed/memory}
#'     \item{"streaming"}{ Processes one trial at a time. Lowest memory usage}
#'   }
#' @param chunk_size Integer number of trials per chunk when using "chunked" strategy
#' @param ram_limit_GB Memory limit in GB for auto strategy selection
#' @param n_cores Number of CPU cores to use for parallel processing (1 = sequential)
#' @param progress Logical whether to show progress bar (requires progressr package)
#' @param verbose Logical whether to print progress messages
#'
#' @return T x V matrix of trial-wise beta estimates
#' @export
run_lss_voxel_loop_core <- function(Y_proj_matrix,
                                   X_trial_onset_list_of_matrices,
                                   B_reconstructor_matrix,
                                   Xi_smoothed_allvox_matrix,
                                   A_lss_fixed_matrix = NULL,
                                   memory_strategy = "auto",
                                   chunk_size = 50,
                                   ram_limit_GB = 4,
                                   n_cores = 1,
                                   progress = TRUE,
                                   verbose = TRUE) {
  
  # Input validation
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  T_trials <- length(X_trial_onset_list_of_matrices)
  p <- ncol(X_trial_onset_list_of_matrices[[1]])
  
  if (nrow(B_reconstructor_matrix) != p) {
    stop("B_reconstructor rows must match trial design matrix columns")
  }
  if (ncol(Xi_smoothed_allvox_matrix) != V) {
    stop("Xi_smoothed columns must match number of voxels")
  }
  
  # Reconstruct HRF shapes
  H_shapes_allvox <- reconstruct_hrf_shapes_core(
    B_reconstructor_matrix,
    Xi_smoothed_allvox_matrix
  )
  
  # Select strategy if auto
  if (memory_strategy == "auto") {
    memory_strategy <- .select_memory_strategy(
      n, V, T_trials, chunk_size, ram_limit_GB, verbose
    )
  }
  
  # Validate strategy
  memory_strategy <- match.arg(
    memory_strategy, 
    c("full", "chunked", "streaming")
  )
  
  if (verbose) {
    message(sprintf(
      "Running LSS with '%s' memory strategy (T=%d trials, V=%d voxels)",
      memory_strategy, T_trials, V
    ))
  }
  
  # Dispatch to appropriate implementation
  result <- switch(memory_strategy,
    full = .lss_execute_full(
      Y_proj_matrix, X_trial_onset_list_of_matrices,
      H_shapes_allvox, A_lss_fixed_matrix, 
      n_cores, progress, verbose
    ),
    chunked = .lss_execute_chunked(
      Y_proj_matrix, X_trial_onset_list_of_matrices,
      H_shapes_allvox, A_lss_fixed_matrix, 
      chunk_size, n_cores, progress, verbose
    ),
    streaming = .lss_execute_streaming(
      Y_proj_matrix, X_trial_onset_list_of_matrices,
      H_shapes_allvox, A_lss_fixed_matrix, 
      n_cores, progress, verbose
    )
  )
  
  return(result)
}

#' Prepare fixed LSS components
#'
#' Precomputes matrices needed for LSS that don't change across trials
#'
#' @param A_fixed_regressors_matrix Optional n x q matrix of fixed nuisance regressors
#' @param lambda_ridge_A Ridge penalty for fixed regressor inversion
#'
#' @return List with P_lss (projection matrix) and has_intercept flag
#' @export
prepare_lss_fixed_components_core <- function(A_fixed_regressors_matrix = NULL,
                                            lambda_ridge_A = 1e-6) {
  if (is.null(A_fixed_regressors_matrix)) {
    return(list(P_lss = NULL, has_intercept = FALSE))
  }
  
  # Check for intercept column (efficient constant detector)
  has_intercept <- any(apply(A_fixed_regressors_matrix, 2, function(x) {
    all(abs(x - x[1]) < 1e-12)
  }))
  
  # fmrilss handles the projection matrix computation internally
  P_lss <- NULL
  
  return(list(
    P_lss = P_lss,
    has_intercept = has_intercept
  ))
}

#' Reconstruct HRF shapes from manifold coordinates
#'
#' @param B_reconstructor_matrix p x m manifold basis matrix
#' @param Xi_manifold_coords_matrix m x V manifold coordinates
#'
#' @return p x V matrix of reconstructed HRF shapes
#' @export
reconstruct_hrf_shapes_core <- function(B_reconstructor_matrix,
                                      Xi_manifold_coords_matrix) {
  # Simple matrix multiplication
  B_reconstructor_matrix %*% Xi_manifold_coords_matrix
}

#' Run LSS for a single voxel
#'
#' Core single-voxel LSS implementation using fmrilss
#'
#' @param Y_proj_voxel_vector Length n vector of data for one voxel
#' @param X_trial_onset_list_of_matrices List of T matrices (n x p each)
#' @param h_voxel_shape_vector Length p HRF shape for this voxel
#' @param A_lss_fixed_matrix Optional n x q fixed regressors
#'
#' @return Length T vector of trial beta estimates
#' @export
run_lss_for_voxel_core <- function(Y_proj_voxel_vector,
                                  X_trial_onset_list_of_matrices,
                                  h_voxel_shape_vector,
                                  A_lss_fixed_matrix = NULL) {
  
  n <- length(Y_proj_voxel_vector)
  T_trials <- length(X_trial_onset_list_of_matrices)
  
  # Create trial regressors by convolving with HRF
  C_matrix <- matrix(0, n, T_trials)
  for (t in seq_len(T_trials)) {
    C_matrix[, t] <- X_trial_onset_list_of_matrices[[t]] %*% h_voxel_shape_vector
  }
  
  # Use fmrilss for the computation
  betas <- fmrilss::lss(
    Y = matrix(Y_proj_voxel_vector, ncol = 1),
    X = C_matrix,
    Z = A_lss_fixed_matrix,
    method = "r_optimized"
  )
  return(as.vector(betas))
}

# ==============================================================================
# Internal Strategy Implementations
# ==============================================================================

#' Execute LSS with full precomputation strategy
#' @noRd
.lss_execute_full <- function(Y_proj_matrix, X_trial_onset_list,
                             H_shapes_allvox, A_lss_fixed_matrix,
                             n_cores, progress, verbose) {
  
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  T_trials <- length(X_trial_onset_list)
  
  if (verbose) {
    message("Precomputing all trial regressors...")
  }
  
  # Precompute all R_t matrices (n x V for each trial)
  R_list <- vector("list", T_trials)
  for (t in seq_len(T_trials)) {
    R_list[[t]] <- X_trial_onset_list[[t]] %*% H_shapes_allvox
  }
  
  if (verbose) message("Running fmrilss on all voxels...")
  
  # Define worker function for a single voxel
  full_worker <- function(v_idx, R_matrices, Y_data, A_fixed) {
    # Extract convolved regressors for this voxel
    C_voxel <- matrix(0, n, T_trials)
    for (t in seq_len(T_trials)) {
      C_voxel[, t] <- R_matrices[[t]][, v_idx]
    }
    
    # Run LSS for this voxel
    betas <- fmrilss::lss(
      Y = Y_data[, v_idx, drop = FALSE],
      X = C_voxel,
      Z = A_fixed,
      method = "r_optimized"
    )
    return(as.vector(betas))
  }
  
  # Process voxels using parallel helper
  voxel_indices <- seq_len(V)
  results_list <- .lss_process_voxels(
    voxel_indices = voxel_indices,
    worker_fun = full_worker,
    R_matrices = R_list,
    Y_data = Y_proj_matrix,
    A_fixed = A_lss_fixed_matrix,
    n_cores = n_cores,
    progress = progress,
    .globals = c("R_list", "Y_proj_matrix", "A_lss_fixed_matrix", "n", "T_trials")
  )
  
  # Combine results into matrix
  Beta_trial_allvox <- do.call(cbind, results_list)
  
  return(Beta_trial_allvox)
}

#' Execute LSS with chunked processing strategy
#' @noRd
.lss_execute_chunked <- function(Y_proj_matrix, X_trial_onset_list,
                                H_shapes_allvox, A_lss_fixed_matrix,
                                chunk_size, n_cores, progress, verbose) {
  
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  T_trials <- length(X_trial_onset_list)
  
  # Calculate chunks
  n_chunks <- ceiling(T_trials / chunk_size)
  Beta_trial_allvox <- matrix(0, T_trials, V)
  
  if (verbose) {
    message(sprintf("Processing %d trials in %d chunks of size %d",
                   T_trials, n_chunks, chunk_size))
  }
  
  # Outer loop over chunks remains sequential
  for (chunk_idx in seq_len(n_chunks)) {
    # Determine trials in this chunk
    trial_start <- (chunk_idx - 1) * chunk_size + 1
    trial_end <- min(chunk_idx * chunk_size, T_trials)
    trial_indices <- trial_start:trial_end
    chunk_trials <- length(trial_indices)
    
    if (verbose) {
      message(sprintf("Processing chunk %d/%d (trials %d-%d)",
                     chunk_idx, n_chunks, trial_start, trial_end))
    }
    
    # Define worker for this chunk
    chunk_worker <- function(v_idx, X_trials_chunk, H_shapes, Y_data, A_fixed, trials_idx) {
      h_voxel <- H_shapes[, v_idx]
      
      # Create design for this voxel and chunk
      C_voxel_chunk <- matrix(0, n, length(trials_idx))
      for (i in seq_along(trials_idx)) {
        t <- trials_idx[i]
        C_voxel_chunk[, i] <- X_trials_chunk[[t]] %*% h_voxel
      }
      
      # Run LSS
      betas <- fmrilss::lss(
        Y = Y_data[, v_idx, drop = FALSE],
        X = C_voxel_chunk,
        Z = A_fixed,
        method = "r_optimized"
      )
      return(as.vector(betas))
    }
    
    # Process voxels in parallel for this chunk
    voxel_indices <- seq_len(V)
    chunk_results <- .lss_process_voxels(
      voxel_indices = voxel_indices,
      worker_fun = chunk_worker,
      X_trials_chunk = X_trial_onset_list,
      H_shapes = H_shapes_allvox,
      Y_data = Y_proj_matrix,
      A_fixed = A_lss_fixed_matrix,
      trials_idx = trial_indices,
      n_cores = n_cores,
      progress = progress && chunk_idx == 1,  # Show progress for first chunk only
      .globals = c("X_trial_onset_list", "H_shapes_allvox", "Y_proj_matrix", 
                   "A_lss_fixed_matrix", "trial_indices", "n")
    )
    
    # Store results for this chunk
    Beta_chunk <- do.call(cbind, chunk_results)
    Beta_trial_allvox[trial_indices, ] <- Beta_chunk
  }
  
  return(Beta_trial_allvox)
}

#' Execute LSS with streaming strategy
#' @noRd
.lss_execute_streaming <- function(Y_proj_matrix, X_trial_onset_list,
                                  H_shapes_allvox, A_lss_fixed_matrix,
                                  n_cores, progress, verbose) {
  
  n <- nrow(Y_proj_matrix)
  V <- ncol(Y_proj_matrix)
  T_trials <- length(X_trial_onset_list)
  
  if (verbose) {
    message("Processing trials one at a time (streaming mode)")
  }
  
  # Define worker for streaming - all computation happens inside
  streaming_worker <- function(v_idx, X_trials, H_shapes, Y_data, A_fixed) {
    # Extract voxel data and HRF
    y_voxel <- Y_data[, v_idx]
    h_voxel <- H_shapes[, v_idx]
    
    # Compute trial regressors on-the-fly for this voxel
    C_voxel <- matrix(0, n, T_trials)
    for (t in seq_len(T_trials)) {
      C_voxel[, t] <- X_trials[[t]] %*% h_voxel
    }
    
    # Run LSS for this voxel
    betas <- fmrilss::lss(
      Y = matrix(y_voxel, ncol = 1),
      X = C_voxel,
      Z = A_fixed,
      method = "r_optimized"
    )
    return(as.vector(betas))
  }
  
  # Process all voxels
  voxel_indices <- seq_len(V)
  results_list <- .lss_process_voxels(
    voxel_indices = voxel_indices,
    worker_fun = streaming_worker,
    X_trials = X_trial_onset_list,
    H_shapes = H_shapes_allvox,
    Y_data = Y_proj_matrix,
    A_fixed = A_lss_fixed_matrix,
    n_cores = n_cores,
    progress = progress,
    .globals = c("X_trial_onset_list", "H_shapes_allvox", "Y_proj_matrix", 
                 "A_lss_fixed_matrix", "n", "T_trials")
  )
  
  # Combine results
  Beta_trial_allvox <- do.call(cbind, results_list)
  
  return(Beta_trial_allvox)
}

# ==============================================================================
# Internal Helper Functions
# ==============================================================================

#' Process voxels in parallel or sequential
#' 
#' Central helper for parallel execution across voxels
#' @param voxel_indices Integer vector of voxel indices to process
#' @param worker_fun Function to apply to each voxel
#' @param ... Additional arguments passed to worker_fun
#' @param n_cores Number of cores to use (1 = sequential)
#' @param progress Show progress bar
#' @param .globals Character vector of global variables to export
#' @param .seed Ensure reproducibility
#' @noRd
.lss_process_voxels <- function(voxel_indices, 
                               worker_fun, 
                               ..., 
                               n_cores = 1, 
                               progress = FALSE,
                               .globals = "auto",
                               .seed = TRUE) {
  
  # Sequential execution (avoid parallel overhead)
  if (n_cores <= 1) {
    if (progress && requireNamespace("progressr", quietly = TRUE)) {
      p <- progressr::progressor(steps = length(voxel_indices))
      results <- lapply(voxel_indices, function(v) {
        res <- worker_fun(v, ...)
        p()
        res
      })
    } else {
      results <- lapply(voxel_indices, worker_fun, ...)
    }
    return(results)
  }
  
  # Parallel execution
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    warning("future.apply not installed, falling back to sequential")
    return(.lss_process_voxels(voxel_indices, worker_fun, ..., 
                              n_cores = 1, progress = progress))
  }
  
  # Set up progress reporting for parallel execution
  if (progress && requireNamespace("progressr", quietly = TRUE)) {
    p <- progressr::progressor(steps = length(voxel_indices))
    worker_fun_with_progress <- function(v, ...) {
      res <- worker_fun(v, ...)
      p()
      res
    }
  } else {
    worker_fun_with_progress <- worker_fun
  }
  
  # Chunk work to reduce overhead
  chunk_size <- ceiling(length(voxel_indices) / n_cores)
  
  # Execute in parallel
  results <- future.apply::future_lapply(
    X = voxel_indices,
    FUN = worker_fun_with_progress,
    ...,
    future.seed = .seed,
    future.globals = .globals,
    future.chunk.size = chunk_size
  )
  
  return(results)
}

#' Select optimal memory strategy based on data size
#' @noRd
.select_memory_strategy <- function(n, V, T_trials, chunk_size, ram_limit_GB, verbose) {
  # Estimate memory requirements
  bytes_per_double <- 8
  
  # Full strategy: stores T matrices of size n x V, plus Y_proj and H_shapes
  # Account for: original data + HRF shapes + convolved regressors + working memory
  full_memory_GB <- ((T_trials + 2) * n * V * bytes_per_double * 1.5) / 1e9
  
  # Chunked strategy: stores chunk_size matrices at a time
  chunked_memory_GB <- ((chunk_size + 2) * n * V * bytes_per_double) / 1e9
  
  # Streaming: minimal memory, processes one trial at a time
  streaming_memory_GB <- (n * V * bytes_per_double) / 1e9
  
  if (verbose) {
    message(sprintf(
      "Memory estimates - Full: %.2f GB, Chunked: %.2f GB, Streaming: %.2f GB",
      full_memory_GB, chunked_memory_GB, streaming_memory_GB
    ))
  }
  
  # Select strategy
  if (full_memory_GB < ram_limit_GB * 0.5) {
    return("full")
  } else if (chunked_memory_GB < ram_limit_GB * 0.7) {
    return("chunked")
  } else {
    return("streaming")
  }
}

# ==============================================================================
# Validation Functions
# ==============================================================================

#' Validate design matrix list
#'
#' @param X_list List of design matrices
#' @param n_timepoints Expected number of timepoints
#' @keywords internal
validate_design_matrix_list <- function(X_list, n_timepoints) {
  if (!is.list(X_list)) {
    stop("X_list must be a list")
  }
  
  if (length(X_list) == 0) {
    stop("X_list cannot be empty")
  }
  
  # Check each matrix
  for (i in seq_along(X_list)) {
    if (!is.matrix(X_list[[i]])) {
      stop(sprintf("Element %d of X_list must be a matrix", i))
    }
    
    if (nrow(X_list[[i]]) != n_timepoints) {
      stop(sprintf(
        "Element %d has %d rows, expected %d",
        i, nrow(X_list[[i]]), n_timepoints
      ))
    }
  }
  
  # Check that all have same number of columns
  n_cols <- vapply(X_list, ncol, integer(1))
  if (length(unique(n_cols)) > 1) {
    stop("All design matrices must have the same number of columns")
  }
  
  invisible(TRUE)
}

#' Validate HRF shape matrix
#'
#' @param H_matrix HRF shape matrix
#' @param n_timepoints Expected number of HRF timepoints
#' @param n_voxels Expected number of voxels
#' @keywords internal
validate_hrf_shape_matrix <- function(H_matrix, n_timepoints, n_voxels) {
  if (!is.matrix(H_matrix)) {
    stop("H_matrix must be a matrix")
  }
  
  if (nrow(H_matrix) != n_timepoints) {
    stop(sprintf(
      "H_matrix has %d rows, expected %d",
      nrow(H_matrix), n_timepoints
    ))
  }
  
  if (ncol(H_matrix) != n_voxels) {
    stop(sprintf(
      "H_matrix has %d columns, expected %d",
      ncol(H_matrix), n_voxels
    ))
  }
  
  invisible(TRUE)
}

# ==============================================================================
# User-facing wrapper functions
# ==============================================================================

#' Run LSS analysis for single voxel (simplified interface)
#'
#' @param y_voxel Timeseries for single voxel
#' @param X_trial_list List of trial onset matrices
#' @param h_voxel HRF shape for the voxel
#' @param TR Repetition time
#' @export
run_lss_for_voxel <- function(y_voxel, X_trial_list, h_voxel, TR = 2) {
  run_lss_for_voxel_core(y_voxel, X_trial_list, h_voxel, NULL)
}

#' Run LSS analysis across all voxels (simplified interface)
#'
#' @param Y_matrix Data matrix (timepoints x voxels)
#' @param X_trial_list List of trial design matrices
#' @param H_matrix HRF shapes (timepoints x voxels)
#' @param memory_strategy Memory optimization strategy
#' @param n_cores Number of CPU cores for parallel processing
#' @param progress Show progress bar
#' @param verbose Print progress messages
#' @export
run_lss_voxel_loop <- function(Y_matrix, X_trial_list, H_matrix,
                              memory_strategy = "auto",
                              n_cores = 1,
                              progress = TRUE,
                              verbose = TRUE) {
  
  # Create dummy manifold components (for compatibility)
  p <- nrow(H_matrix)
  m <- min(5, p)  # Use small manifold dimension
  
  # Create identity-like reconstructor
  B_reconstructor <- matrix(0, p, m)
  for (i in seq_len(min(m, p))) {
    B_reconstructor[i, i] <- 1
  }
  
  # Project HRFs to get coordinates (m x V matrix)
  Xi_coords <- solve(crossprod(B_reconstructor) + 1e-6 * diag(m)) %*% 
               t(B_reconstructor) %*% H_matrix
  
  # Validate dimensions
  stopifnot(nrow(Xi_coords) == m, ncol(Xi_coords) == ncol(H_matrix))
  
  # Call main function
  run_lss_voxel_loop_core(
    Y_matrix, X_trial_list, B_reconstructor, Xi_coords,
    A_lss_fixed_matrix = NULL,
    memory_strategy = memory_strategy,
    n_cores = n_cores,
    progress = progress,
    verbose = verbose
  )
}

#' Fast LSS implementation
#'
#' Wrapper for fmrilss package functionality
#'
#' @param Y Data matrix
#' @param dmat_base Base design matrix
#' @param dmat_ran Random effects design
#' @param dmat_fixed Fixed effects design
#' @export
lss_fast <- function(Y, dmat_base, dmat_ran, dmat_fixed = NULL) {
  # Combine base and random designs
  X_combined <- cbind(dmat_base, dmat_ran)
  
  # Call fmrilss
  fmrilss::lss(
    Y = Y,
    X = X_combined,
    Z = dmat_fixed,
    method = "r_optimized"
  )
}