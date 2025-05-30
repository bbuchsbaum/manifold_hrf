# Core Manifold Construction Functions (Component 0)
# Implementation of MHRF-CORE-MANIFOLD-01 and MHRF-CORE-MANIFOLD-02

#' Calculate Manifold Affinity and Markov Matrix (Core)
#'
#' Computes a Markov transition matrix for HRF manifold construction using
#' self-tuning bandwidth for local scaling of affinities.
#'
#' @param L_library_matrix A p x N matrix of HRF shapes, where p is the number 
#'   of time points and N is the number of HRFs in the library
#' @param k_local_nn_for_sigma Integer, k-nearest neighbors for self-tuning 
#'   bandwidth calculation (e.g., 7)
#' @param use_sparse_W_params List with optional parameters for sparse W matrix:
#'   \itemize{
#'     \item \code{sparse_if_N_gt}: Threshold for N to switch to sparse matrix (e.g., 5000)
#'     \item \code{k_nn_for_W_sparse}: Number of nearest neighbors to keep in sparse W
#'   }
#' 
#' @return S_markov_matrix An N x N Markov transition matrix (regular or sparse 
#'   Matrix format depending on parameters). Each row sums to 1.
#'   
#' @details This function implements the affinity matrix construction from
#'   Component 0, Step 1 of the M-HRF-LSS pipeline. It uses the self-tuning
#'   local scaling method from Zelnik-Manor & Perona (2005) where each HRF's
#'   bandwidth sigma_i is set to its k-th nearest neighbor distance.
#'   
#' @examples
#' \dontrun{
#' # Create synthetic HRF library
#' p <- 30  # time points
#' N <- 100 # number of HRFs
#' L_library <- matrix(rnorm(p * N), nrow = p, ncol = N)
#' 
#' # Compute Markov matrix
#' S <- calculate_manifold_affinity_core(L_library, k_local_nn_for_sigma = 7)
#' 
#' # With sparse matrix for large N
#' S_sparse <- calculate_manifold_affinity_core(
#'   L_library, 
#'   k_local_nn_for_sigma = 7,
#'   use_sparse_W_params = list(sparse_if_N_gt = 50, k_nn_for_W_sparse = 20)
#' )
#' }
#' 
#' @export
calculate_manifold_affinity_core <- function(L_library_matrix,
                                            k_local_nn_for_sigma,
                                            use_sparse_W_params = list()) {
  
  # Input validation
  if (!is.matrix(L_library_matrix)) {
    stop("L_library_matrix must be a matrix")
  }
  
  N <- ncol(L_library_matrix)
  p <- nrow(L_library_matrix)
  
  if (k_local_nn_for_sigma >= N) {
    stop("k_local_nn_for_sigma must be less than the number of HRFs (N)")
  }
  
  # Extract parameters for sparse matrix handling
  sparse_threshold <- use_sparse_W_params$sparse_if_N_gt
  k_nn_sparse <- use_sparse_W_params$k_nn_for_W_sparse
  
  # Step 1: Compute pairwise Euclidean distances
  # Using dist() for efficiency, then convert to matrix
  dist_mat <- as.matrix(dist(t(L_library_matrix)))
  
  # Step 2: Find k-th nearest neighbor distance for each HRF (self-tuning bandwidth)
  # For each HRF, sort distances and take the k-th smallest (excluding self)
  sigma_i <- numeric(N)
  for (i in 1:N) {
    # Sort distances from HRF i (excluding self-distance which is 0)
    sorted_dists <- sort(dist_mat[i, -i])
    # k-th nearest neighbor distance
    sigma_i[i] <- sorted_dists[k_local_nn_for_sigma]
    
    # Ensure sigma is not zero (can happen with duplicate HRFs)
    if (sigma_i[i] == 0) {
      # Use median of all non-zero distances as fallback
      non_zero_dists <- sorted_dists[sorted_dists > 0]
      if (length(non_zero_dists) > 0) {
        sigma_i[i] <- median(non_zero_dists)
      } else {
        # Last resort: use a small positive value
        sigma_i[i] <- 1e-6
      }
    }
  }
  
  # Step 3: Compute affinity matrix W with self-tuning bandwidth
  # W_ij = exp(-dist_ij^2 / (sigma_i * sigma_j))
  W <- matrix(0, N, N)
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        W[i, j] <- exp(-dist_mat[i, j]^2 / (sigma_i[i] * sigma_i[j]))
      }
    }
  }
  # Set diagonal to 0 (no self-loops)
  diag(W) <- 0
  
  # Step 4: Optional sparsification for large N
  if (!is.null(sparse_threshold) && N > sparse_threshold && !is.null(k_nn_sparse)) {
    # Convert to sparse matrix keeping only k_nn_sparse nearest neighbors
    W_sparse <- W
    for (i in 1:N) {
      # Find indices of k_nn_sparse largest values (nearest neighbors in affinity space)
      top_k_indices <- order(W[i, ], decreasing = TRUE)[1:k_nn_sparse]
      # Zero out all other entries
      mask <- rep(FALSE, N)
      mask[top_k_indices] <- TRUE
      W_sparse[i, !mask] <- 0
    }
    # Symmetrize the sparse matrix (take maximum of W_ij and W_ji)
    W <- pmax(W_sparse, t(W_sparse))
    
    # Convert to sparse matrix format
    if (requireNamespace("Matrix", quietly = TRUE)) {
      W <- Matrix::Matrix(W, sparse = TRUE)
    }
  }
  
  # Step 5: Create Markov matrix S = D^(-1) * W
  # D is the degree matrix (diagonal matrix of row sums)
  if (inherits(W, "Matrix")) {
    row_sums <- Matrix::rowSums(W)
  } else {
    row_sums <- rowSums(W)
  }
  
  # Handle potential zero row sums (isolated nodes) or NA values
  if (any(is.na(row_sums)) || any(row_sums == 0)) {
    if (any(is.na(row_sums))) {
      warning("Some row sums are NA. Setting to 1.")
      row_sums[is.na(row_sums)] <- 1
    }
    if (any(row_sums == 0)) {
      warning("Some HRFs have zero affinity to all others. Setting their row to uniform distribution.")
      row_sums[row_sums == 0] <- 1
    }
  }
  
  # Create D_inv (inverse degree matrix)
  D_inv <- diag(1 / row_sums)
  
  # Compute Markov matrix S
  if (inherits(W, "Matrix")) {
    # If W is sparse, keep S sparse
    D_inv <- Matrix::Diagonal(x = 1 / row_sums)
    S_markov_matrix <- D_inv %*% W
  } else {
    S_markov_matrix <- D_inv %*% W
  }
  
  return(S_markov_matrix)
}

#' Get Manifold Basis Reconstructor (Core)
#'
#' Computes diffusion map coordinates and HRF reconstructor matrix for the manifold.
#'
#' @param S_markov_matrix N x N Markov matrix from calculate_manifold_affinity_core
#' @param L_library_matrix p x N matrix of HRF shapes (same as used for S_markov_matrix)
#' @param m_manifold_dim_target Target manifold dimensionality (e.g., 3-5)
#' @param m_manifold_dim_min_variance Minimum variance explained threshold (default 0.95)
#' 
#' @return A list containing:
#'   \itemize{
#'     \item \code{B_reconstructor_matrix}: p x m matrix mapping manifold coords to HRF shapes
#'     \item \code{Phi_coords_matrix}: N x m matrix of diffusion map coordinates
#'     \item \code{eigenvalues_S_vector}: Vector of eigenvalues from S decomposition
#'     \item \code{m_final_dim}: Final manifold dimension used
#'     \item \code{m_auto_selected_dim}: Automatically selected dimension based on variance
#'   }
#'   
#' @details This function implements Component 0, Steps 3-5 of the M-HRF-LSS pipeline.
#'   It computes the diffusion map embedding of the HRF library and creates a 
#'   reconstructor matrix that maps low-dimensional manifold coordinates back to
#'   full HRF shapes.
#'   
#' @examples
#' \dontrun{
#' # Create synthetic HRF library and compute manifold
#' p <- 30
#' N <- 100
#' L_library <- matrix(rnorm(p * N), nrow = p, ncol = N)
#' S <- calculate_manifold_affinity_core(L_library, k_local_nn_for_sigma = 7)
#' 
#' # Get manifold basis
#' manifold <- get_manifold_basis_reconstructor_core(
#'   S, L_library, 
#'   m_manifold_dim_target = 5,
#'   m_manifold_dim_min_variance = 0.95
#' )
#' }
#' 
#' @export
get_manifold_basis_reconstructor_core <- function(S_markov_matrix,
                                                 L_library_matrix,
                                                 m_manifold_dim_target,
                                                 m_manifold_dim_min_variance = 0.95) {
  if (!is.matrix(S_markov_matrix) && !inherits(S_markov_matrix, "Matrix")) {
    stop("S_markov_matrix must be a matrix or Matrix object")
  }
  
  if (!is.matrix(L_library_matrix)) {
    stop("L_library_matrix must be a matrix")
  }
  
  N <- ncol(L_library_matrix)
  p <- nrow(L_library_matrix)
  
  if (nrow(S_markov_matrix) != N || ncol(S_markov_matrix) != N) {
    stop("S_markov_matrix dimensions must match number of HRFs in L_library_matrix")
  }
  
  if (m_manifold_dim_min_variance <= 0 || m_manifold_dim_min_variance > 1) {
    stop("m_manifold_dim_min_variance must be between 0 and 1")
  }
  
  # Step 3: Compute diffusion map coordinates using eigendecomposition
  # We need top k eigenvectors, where k is larger than target to allow for selection
  k_eig_max <- min(N - 1, max(10, m_manifold_dim_target + 5))
  
  # Use RSpectra for efficient eigendecomposition of large matrices
  if (requireNamespace("RSpectra", quietly = TRUE) && N > 100) {
    # For large matrices, use RSpectra
    eig_result <- RSpectra::eigs_sym(
      S_markov_matrix, 
      k = k_eig_max + 1,  # +1 for the trivial eigenvector
      which = "LM"  # Largest magnitude
    )
    eigenvalues_full <- eig_result$values
    eigenvectors_full <- eig_result$vectors
  } else {
    # For smaller matrices, use base R eigen
    if (inherits(S_markov_matrix, "Matrix")) {
      S_markov_matrix <- as.matrix(S_markov_matrix)
    }
    eig_result <- eigen(S_markov_matrix, symmetric = FALSE)
    # Take only the top k_eig_max + 1 eigenvectors
    eigenvalues_full <- eig_result$values[1:(k_eig_max + 1)]
    eigenvectors_full <- eig_result$vectors[, 1:(k_eig_max + 1)]
  }
  
  # The first eigenvector should be trivial (all ones for stochastic matrix)
  # We exclude it from the manifold coordinates
  eigenvalues_S <- eigenvalues_full[-1]
  Phi_raw_full <- eigenvectors_full[, -1, drop = FALSE]
  
  # Automatic m selection based on variance explained
  # Using absolute values of eigenvalues for variance calculation
  abs_eigenvalues <- abs(eigenvalues_S)
  cum_var_explained <- cumsum(abs_eigenvalues) / sum(abs_eigenvalues)
  m_auto <- which(cum_var_explained >= m_manifold_dim_min_variance)[1]
  
  # If no dimension meets the variance threshold, use all available
  if (is.na(m_auto)) {
    m_auto <- length(eigenvalues_S)
    warning(sprintf("Could not achieve %.1f%% variance with available dimensions. Using all %d dimensions.", 
                   m_manifold_dim_min_variance * 100, m_auto))
  }
  
  # Choose final dimension (use target dimension, but limit to available eigenvectors)
  m_final <- min(m_manifold_dim_target, length(eigenvalues_S))

  # Warn if target is much lower than auto-selected dimension needed for variance threshold
  if (m_manifold_dim_target < m_auto * 0.5) {
    warning(sprintf("Target dimension %d is much lower than auto-selected %d (for %.1f%% variance)", 
                   m_manifold_dim_target, m_auto, m_manifold_dim_min_variance * 100))
  }
  
  # Step 4: Extract final manifold coordinates
  Phi_coords_matrix <- Phi_raw_full[, 1:m_final, drop = FALSE]
  
  # Enforce consistent sign for reproducibility 
  # (first non-zero element of each eigenvector should be positive)
  for (j in 1:m_final) {
    first_nonzero_idx <- which(abs(Phi_coords_matrix[, j]) > 1e-10)[1]
    if (!is.na(first_nonzero_idx) && Phi_coords_matrix[first_nonzero_idx, j] < 0) {
      Phi_coords_matrix[, j] <- -Phi_coords_matrix[, j]
    }
  }
  
  # Step 5: Compute HRF reconstructor matrix B
  # B = L * Phi * (Phi' * Phi + ridge * I)^(-1)
  # This gives us the best linear reconstruction of HRFs from manifold coordinates
  
  # Ridge parameter for numerical stability
  ridge_param <- 1e-8
  
  # Compute Phi' * Phi
  PhiTPhi <- crossprod(Phi_coords_matrix)
  
  # Add ridge regularization
  PhiTPhi_reg <- PhiTPhi + ridge_param * diag(m_final)
  
  # Compute reconstructor
  B_reconstructor_matrix <- L_library_matrix %*% Phi_coords_matrix %*% solve(PhiTPhi_reg)
  
  # Return results
  list(
    B_reconstructor_matrix = B_reconstructor_matrix,
    Phi_coords_matrix = Phi_coords_matrix,
    eigenvalues_S_vector = eigenvalues_full,  # Return all eigenvalues for diagnostics
    m_final_dim = m_final,
    m_auto_selected_dim = m_auto,
    m_manifold_dim = m_final
  )
}
