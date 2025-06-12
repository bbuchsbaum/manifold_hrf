# Core Manifold Construction Functions (Component 0) - REFACTORED
# Addresses critical issues from code review

#' Build Sparse Affinity Matrix from HRF Library
#'
#' Computes a sparse affinity matrix using self-tuning local scaling,
#' with proper handling of normalization and isolated nodes.
#'
#' @param L_library_matrix A p x N matrix of HRF shapes, where p is the number 
#'   of time points and N is the number of HRFs in the library
#' @param k Number of nearest neighbors for graph construction
#' @param k_sigma Number of nearest neighbors for sigma estimation (default = k)
#' @param normalization Normalization method: "symmetric" (for diffusion maps) 
#'   or "row_stochastic" (for Markov chains). Default is "symmetric".
#' @param handle_isolated How to handle isolated nodes: "error", "warn_remove", 
#'   or "connect_to_self". Default is "warn_remove".
#' @param ann_method Backend for k-NN search: "auto", "RANN", "RcppHNSW", or "exact"
#' @param return_diagnostics Whether to return additional diagnostic information
#' 
#' @return A list object of class 'manifold_affinity' containing:
#'   \itemize{
#'     \item \code{W}: Raw sparse affinity matrix
#'     \item \code{T_norm}: Normalized sparse matrix (symmetric or row-stochastic)
#'     \item \code{normalization}: Normalization method used
#'     \item \code{params}: List of parameters used
#'     \item \code{diagnostics}: Optional diagnostics (sigma values, isolated nodes, etc.)
#'   }
#'   
#' @details This function implements a corrected version of the affinity matrix 
#'   construction from Component 0 of the M-HRF-LSS pipeline. Key improvements:
#'   - Always uses sparse matrices to avoid memory issues
#'   - Properly handles normalization (symmetric vs row-stochastic)
#'   - No global diagonal regularization (only isolated node handling)
#'   - Returns rich diagnostic information
#'   
#' @export
build_manifold_affinity <- function(L_library_matrix,
                                   k,
                                   k_sigma = NULL,
                                   normalization = c("symmetric", "row_stochastic"),
                                   handle_isolated = c("warn_remove", "error", "connect_to_self"),
                                   ann_method = c("auto", "RANN", "RcppHNSW", "exact"),
                                   return_diagnostics = TRUE) {
  
  # Input validation
  if (!is.matrix(L_library_matrix)) {
    stop("L_library_matrix must be a matrix")
  }
  
  N <- ncol(L_library_matrix)
  p <- nrow(L_library_matrix)
  
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k != round(k)) {
    stop("k must be a positive integer")
  }
  
  if (k >= N) {
    stop("k must be less than the number of HRFs (N)")
  }
  
  # Default k_sigma to k if not specified
  if (is.null(k_sigma)) {
    k_sigma <- k
  } else {
    if (!is.numeric(k_sigma) || length(k_sigma) != 1 || 
        k_sigma < 1 || k_sigma != round(k_sigma)) {
      stop("k_sigma must be a positive integer")
    }
    if (k_sigma >= N) {
      stop("k_sigma must be less than the number of HRFs (N)")
    }
  }
  
  normalization <- match.arg(normalization)
  handle_isolated <- match.arg(handle_isolated)
  ann_method <- match.arg(ann_method)
  
  # Timing for diagnostics
  timing <- list()
  t0 <- Sys.time()
  
  # Step 1: k-NN search (always sparse approach)
  # Choose backend based on availability and preference
  knn_backend <- NULL
  
  if (ann_method == "auto") {
    # Auto-select based on availability and N
    if (N > 5000 && requireNamespace("RcppHNSW", quietly = TRUE)) {
      ann_method <- "RcppHNSW"
    } else if (requireNamespace("RANN", quietly = TRUE)) {
      ann_method <- "RANN"
    } else if (exists("knn_search_cpp", mode = "function")) {
      ann_method <- "exact"
    } else {
      stop("No k-NN backend available. Install RANN or RcppHNSW package.")
    }
  }
  
  # Perform k-NN search
  k_search <- max(k, k_sigma)  # Search for max needed
  
  if (ann_method == "RcppHNSW") {
    if (!requireNamespace("RcppHNSW", quietly = TRUE)) {
      stop("RcppHNSW package not available")
    }
    knn_backend <- "RcppHNSW"
    ann_res <- RcppHNSW::hnsw_knn(t(L_library_matrix), k = k_search + 1)
    nn_idx <- ann_res$idx[, -1, drop = FALSE]  # Remove self
    nn_dist <- ann_res$dist[, -1, drop = FALSE]
    
  } else if (ann_method == "RANN") {
    if (!requireNamespace("RANN", quietly = TRUE)) {
      stop("RANN package not available")
    }
    knn_backend <- "RANN"
    nn_res <- RANN::nn2(t(L_library_matrix), k = k_search + 1)
    nn_idx <- nn_res$nn.idx[, -1, drop = FALSE]
    nn_dist <- nn_res$nn.dists[, -1, drop = FALSE]
    
  } else if (ann_method == "exact") {
    if (exists("knn_search_cpp", mode = "function")) {
      knn_backend <- "knn_search_cpp"
      nn_res <- knn_search_cpp(L_library_matrix, L_library_matrix, k_search + 1)
      nn_idx <- t(nn_res$idx)[, -1, drop = FALSE]
      nn_dist <- t(nn_res$dist)[, -1, drop = FALSE]
    } else {
      stop("knn_search_cpp function not found. Ensure package is properly compiled.")
    }
  }
  
  timing$knn <- as.numeric(Sys.time() - t0, units = "secs")
  t0 <- Sys.time()
  
  # Step 2: Extract sigma values (distance to k_sigma-th neighbor)
  sigma_i <- nn_dist[, k_sigma]
  
  # Handle zero sigmas (identical points)
  zero_sigma <- which(sigma_i == 0)
  if (length(zero_sigma) > 0) {
    warning(sprintf("Found %d points with zero k-th neighbor distance. Using median distance.", 
                   length(zero_sigma)))
    # Use median of non-zero distances for these points
    for (i in zero_sigma) {
      nz_dists <- nn_dist[i, nn_dist[i, ] > 0]
      sigma_i[i] <- if (length(nz_dists) > 0) median(nz_dists) else 1e-6
    }
  }
  
  # Step 3: Build sparse affinity matrix W
  # Only store k nearest neighbors per point
  i_vec <- rep(seq_len(N), each = k)
  j_vec <- as.vector(t(nn_idx[, 1:k]))
  d_vec <- as.vector(t(nn_dist[, 1:k]))
  
  # Self-tuning Gaussian kernel
  sigma_prod_vec <- sigma_i[i_vec] * sigma_i[j_vec]
  w_vec <- exp(-(d_vec^2) / sigma_prod_vec)
  
  # Create sparse matrix
  W <- Matrix::sparseMatrix(
    i = i_vec, 
    j = j_vec, 
    x = w_vec,
    dims = c(N, N),
    repr = "T"  # Triplet form for efficiency
  )
  
  # Symmetrize by taking maximum
  W <- Matrix::forceSymmetric(pmax(W, Matrix::t(W)), uplo = "U")
  
  timing$build_W <- as.numeric(Sys.time() - t0, units = "secs")
  t0 <- Sys.time()
  
  # Step 4: Handle isolated nodes
  degrees <- Matrix::rowSums(W)
  isolated_idx <- which(degrees == 0)
  n_isolated <- length(isolated_idx)
  
  if (n_isolated > 0) {
    msg <- sprintf("Found %d isolated nodes (%.1f%% of total)", 
                   n_isolated, 100 * n_isolated / N)
    
    if (handle_isolated == "error") {
      stop(msg)
    } else if (handle_isolated == "warn_remove") {
      warning(msg, ". These will be removed from the manifold.")
      # Mark for removal later
    } else if (handle_isolated == "connect_to_self") {
      warning(msg, ". Connecting to self (creates absorbing states).")
      # Add self-loops
      W <- W + Matrix::Diagonal(N, x = ifelse(seq_len(N) %in% isolated_idx, 1, 0))
      degrees <- Matrix::rowSums(W)  # Recompute
    }
  }
  
  # Step 5: Normalize the affinity matrix
  if (handle_isolated != "warn_remove" || n_isolated == 0) {
    # Normal case: normalize all nodes
    if (normalization == "symmetric") {
      # D^(-1/2) * W * D^(-1/2) for diffusion maps
      # Add small epsilon for numerical stability in sqrt
      D_inv_sqrt <- Matrix::Diagonal(x = 1 / sqrt(degrees + .Machine$double.eps))
      T_norm <- D_inv_sqrt %*% W %*% D_inv_sqrt
      # Ensure symmetry
      T_norm <- Matrix::forceSymmetric(T_norm, uplo = "U")
      
    } else if (normalization == "row_stochastic") {
      # D^(-1) * W for random walk
      D_inv <- Matrix::Diagonal(x = 1 / (degrees + .Machine$double.eps))
      T_norm <- D_inv %*% W
    }
  } else {
    # Remove isolated nodes before normalization
    keep_idx <- setdiff(seq_len(N), isolated_idx)
    W_clean <- W[keep_idx, keep_idx]
    degrees_clean <- degrees[keep_idx]
    
    if (normalization == "symmetric") {
      D_inv_sqrt <- Matrix::Diagonal(x = 1 / sqrt(degrees_clean + .Machine$double.eps))
      T_norm <- D_inv_sqrt %*% W_clean %*% D_inv_sqrt
      T_norm <- Matrix::forceSymmetric(T_norm, uplo = "U")
    } else {
      D_inv <- Matrix::Diagonal(x = 1 / (degrees_clean + .Machine$double.eps))
      T_norm <- D_inv %*% W_clean
    }
  }
  
  timing$normalize <- as.numeric(Sys.time() - t0, units = "secs")
  
  # Prepare result
  result <- structure(
    list(
      W = W,
      T_norm = T_norm,
      normalization = normalization,
      params = list(
        k = k,
        k_sigma = k_sigma,
        handle_isolated = handle_isolated,
        N = N,
        p = p
      )
    ),
    class = "manifold_affinity"
  )
  
  # Add diagnostics if requested
  if (return_diagnostics) {
    result$diagnostics <- list(
      n_isolated = n_isolated,
      isolated_indices = isolated_idx,
      sigma_values = sigma_i,
      knn_backend = knn_backend,
      sparsity = 1 - Matrix::nnzero(W) / (N * N),
      timing_sec = timing
    )
  }
  
  return(result)
}

#' Compute Diffusion Map Basis from Affinity Matrix
#'
#' Computes the diffusion map embedding from a manifold affinity object,
#' with proper handling of the symmetric eigendecomposition.
#'
#' @param affinity An object of class 'manifold_affinity' from build_manifold_affinity
#' @param n_dims Number of dimensions to compute. Can be an integer or "auto"
#' @param min_variance For n_dims="auto", minimum cumulative variance to retain
#' @param return_basis Whether to compute the HRF reconstructor matrix
#' @param L_library_matrix Required if return_basis=TRUE. Original HRF library matrix
#' 
#' @return A list object of class 'diffusion_basis' containing:
#'   \itemize{
#'     \item \code{eigenvalues}: Vector of eigenvalues
#'     \item \code{eigenvectors}: Matrix of eigenvectors (N x n_dims)
#'     \item \code{n_dims}: Number of dimensions used
#'     \item \code{B_reconstructor}: Optional HRF reconstructor matrix (p x n_dims)
#'     \item \code{affinity}: Original affinity object for provenance
#'   }
#'   
#' @details This function correctly handles the eigendecomposition based on
#'   the normalization method used in the affinity matrix:
#'   - For symmetric normalization: Uses RSpectra::eigs_sym for efficiency
#'   - For row-stochastic: Uses standard eigendecomposition
#'   - Always checks for complex eigenvalues and warns appropriately
#'   
#' @export
compute_diffusion_basis <- function(affinity,
                                   n_dims = "auto",
                                   min_variance = 0.95,
                                   return_basis = TRUE,
                                   L_library_matrix = NULL) {
  
  # Validate input
  if (!inherits(affinity, "manifold_affinity")) {
    stop("affinity must be an object from build_manifold_affinity()")
  }
  
  # Extract components
  T_norm <- affinity$T_norm
  normalization <- affinity$normalization
  N <- nrow(T_norm)
  
  # Handle removed isolated nodes
  if (ncol(T_norm) < affinity$params$N) {
    # Some nodes were removed
    N_original <- affinity$params$N
    removed_nodes <- affinity$diagnostics$isolated_indices
    message(sprintf("Working with %d nodes (%d isolated nodes removed)", 
                   N, length(removed_nodes)))
  } else {
    N_original <- N
    removed_nodes <- integer(0)
  }
  
  # Determine number of eigenvectors to compute
  if (is.character(n_dims) && n_dims == "auto") {
    # Compute enough to determine dimensionality
    k_compute <- min(N - 1, max(20, ceiling(N * 0.1)))
  } else {
    if (!is.numeric(n_dims) || length(n_dims) != 1 || 
        n_dims < 1 || n_dims != round(n_dims)) {
      stop("n_dims must be a positive integer or 'auto'")
    }
    k_compute <- min(n_dims + 5, N - 1)  # Compute a few extra
  }
  
  # Eigendecomposition
  timing_eig_start <- Sys.time()
  
  if (N > 100 && requireNamespace("RSpectra", quietly = TRUE)) {
    # Use iterative solver for large matrices
    if (normalization == "symmetric") {
      # Symmetric case - guaranteed real eigenvalues
      eig_result <- RSpectra::eigs_sym(
        T_norm,
        k = k_compute,
        which = "LA"  # Largest algebraic (for non-negative matrix)
      )
      eigenvalues <- eig_result$values
      eigenvectors <- eig_result$vectors
      
    } else {
      # Row-stochastic case - may have complex eigenvalues
      eig_result <- RSpectra::eigs(
        T_norm,
        k = k_compute,
        which = "LM"  # Largest magnitude
      )
      
      # Check for complex eigenvalues
      if (any(abs(Im(eig_result$values)) > 1e-6)) {
        warning("Complex eigenvalues detected. Taking real parts.")
      }
      eigenvalues <- Re(eig_result$values)
      eigenvectors <- Re(eig_result$vectors)
    }
  } else {
    # Small matrices - use base eigen
    T_dense <- as.matrix(T_norm)
    eig_result <- eigen(T_dense, symmetric = (normalization == "symmetric"))
    
    if (normalization != "symmetric" && any(abs(Im(eig_result$values)) > 1e-6)) {
      warning("Complex eigenvalues detected. Taking real parts.")
    }
    
    eigenvalues <- Re(eig_result$values[1:k_compute])
    eigenvectors <- Re(eig_result$vectors[, 1:k_compute, drop = FALSE])
  }
  
  timing_eig <- as.numeric(Sys.time() - timing_eig_start, units = "secs")
  
  # Remove trivial eigenvector (should be first for row-stochastic)
  # For symmetric normalization, there may not be a clear trivial eigenvector
  if (normalization == "row_stochastic") {
    # First eigenvector should be constant
    const_score <- var(eigenvectors[, 1])
    if (const_score < 1e-6) {
      # Remove trivial
      eigenvalues <- eigenvalues[-1]
      eigenvectors <- eigenvectors[, -1, drop = FALSE]
    } else {
      warning("First eigenvector is not constant. Markov matrix may be malformed.")
    }
  }
  
  # Determine final dimensionality
  if (is.character(n_dims) && n_dims == "auto") {
    dim_info <- select_manifold_dim(eigenvalues, min_variance)
    n_dims_final <- dim_info$m_auto
  } else {
    n_dims_final <- min(n_dims, ncol(eigenvectors))
  }
  
  # Extract final eigenvectors
  eigenvectors_final <- eigenvectors[, 1:n_dims_final, drop = FALSE]
  eigenvalues_final <- eigenvalues[1:n_dims_final]
  
  # Enforce consistent sign
  for (j in 1:n_dims_final) {
    first_nonzero <- which(abs(eigenvectors_final[, j]) > 1e-10)[1]
    if (!is.na(first_nonzero) && eigenvectors_final[first_nonzero, j] < 0) {
      eigenvectors_final[, j] <- -eigenvectors_final[, j]
    }
  }
  
  # Compute reconstructor if requested
  B_reconstructor <- NULL
  if (return_basis) {
    if (is.null(L_library_matrix)) {
      stop("L_library_matrix required when return_basis=TRUE")
    }
    
    # Handle removed nodes
    if (length(removed_nodes) > 0) {
      # Use only the kept nodes
      keep_idx <- setdiff(seq_len(N_original), removed_nodes)
      L_kept <- L_library_matrix[, keep_idx, drop = FALSE]
    } else {
      L_kept <- L_library_matrix
    }
    
    # B = L * Phi * (Phi' * Phi + ridge * I)^(-1)
    ridge_param <- 1e-8
    PhiTPhi <- crossprod(eigenvectors_final)
    PhiTPhi_reg <- PhiTPhi + ridge_param * diag(n_dims_final)
    
    # Use Cholesky for SPD system
    chol_PhiTPhi <- chol(PhiTPhi_reg)
    B_reconstructor <- L_kept %*% eigenvectors_final %*% 
                       backsolve(chol_PhiTPhi, backsolve(chol_PhiTPhi, diag(n_dims_final), 
                                                         transpose = TRUE))
  }
  
  # Prepare result
  result <- structure(
    list(
      eigenvalues = eigenvalues_final,
      eigenvectors = eigenvectors_final,
      n_dims = n_dims_final,
      B_reconstructor = B_reconstructor,
      affinity = affinity,
      removed_nodes = removed_nodes,
      timing_eig_sec = timing_eig
    ),
    class = "diffusion_basis"
  )
  
  return(result)
}

# Print methods for the new classes
#' @export
print.manifold_affinity <- function(x, ...) {
  cat("Manifold Affinity Matrix\n")
  cat(sprintf("  Nodes: %d\n", nrow(x$T_norm)))
  cat(sprintf("  Normalization: %s\n", x$normalization))
  cat(sprintf("  k-NN: %d (sigma: %d)\n", x$params$k, x$params$k_sigma))
  if (!is.null(x$diagnostics)) {
    cat(sprintf("  Isolated nodes: %d\n", x$diagnostics$n_isolated))
    cat(sprintf("  Sparsity: %.1f%%\n", 100 * x$diagnostics$sparsity))
    cat(sprintf("  Backend: %s\n", x$diagnostics$knn_backend))
  }
  invisible(x)
}

#' @export
print.diffusion_basis <- function(x, ...) {
  cat("Diffusion Map Basis\n")
  cat(sprintf("  Dimensions: %d\n", x$n_dims))
  cat(sprintf("  Nodes: %d\n", nrow(x$eigenvectors)))
  if (length(x$removed_nodes) > 0) {
    cat(sprintf("  Removed nodes: %d\n", length(x$removed_nodes)))
  }
  cat(sprintf("  Eigenvalues: %s\n", 
              paste(sprintf("%.3f", x$eigenvalues[1:min(5, x$n_dims)]), collapse = " ")))
  if (!is.null(x$B_reconstructor)) {
    cat(sprintf("  Reconstructor: %d x %d\n", 
                nrow(x$B_reconstructor), ncol(x$B_reconstructor)))
  }
  invisible(x)
}