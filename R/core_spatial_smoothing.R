# Core Spatial Smoothing Functions (Component 2)
# Implementation of MHRF-CORE-SPSMOOTH-01 and MHRF-CORE-SPSMOOTH-02

#' Construct Graph Laplacian for Voxel Coordinates with Weights (Internal)
#'
#' Creates a sparse graph Laplacian matrix for spatial smoothing based on
#' k-nearest neighbors in 3D voxel space. This is an internal version that
#' supports weighted graphs.
#'
#' @param voxel_coords_matrix A V x 3 matrix of voxel coordinates, where V is 
#'   the number of voxels and columns are x, y, z coordinates
#' @param num_neighbors_Lsp Number of nearest neighbors for graph construction
#'   (e.g., 6 for face-connected, 18 for edge-connected, 26 for corner-connected)
#' @param distance_engine Choice of distance computation engine. "euclidean"
#'   performs an exact search while "ann_euclidean" uses approximate nearest
#'   neighbors via the \pkg{RcppHNSW} package. If \pkg{RcppHNSW} is not available,
#'   the function falls back to the exact search and issues a warning.
#' @param ann_threshold When \code{distance_engine = "euclidean"}, use the
#'   approximate search for datasets larger than this threshold if
#'   \pkg{RcppHNSW} is installed.
#' @param weight_scheme Weight scheme for the adjacency matrix. "binary" uses
#'   0/1 weights, "gaussian" uses exp(-d²/σ²) where σ is the median distance
#'   to the k-th nearest neighbor across all voxels.
#'   
#' @return L_sp_sparse_matrix A V x V sparse Laplacian matrix (Matrix::sparseMatrix)
#'   
#' @details This function implements Component 2, Step 1 of the M-HRF-LSS pipeline.
#'   It constructs a k-nearest neighbor graph from voxel coordinates and computes
#'   the combinatorial graph Laplacian L = D - W, where W is the symmetrized 
#'   binary adjacency matrix and D is the degree matrix. The initial k-NN graph
#'   is directed (each voxel points to its k nearest neighbors), which is then
#'   symmetrized by taking the maximum of W and its transpose. The resulting 
#'   Laplacian is used for spatial smoothing of manifold coordinates.
#'
#' @examples
#' \dontrun{
#' # Create example 3D grid of voxels
#' coords <- expand.grid(x = 1:10, y = 1:10, z = 1:5)
#' voxel_coords <- as.matrix(coords)
#' 
#' # Construct 6-neighbor Laplacian (face-connected)
#' L_sparse <- .make_voxel_graph_laplacian_weighted(voxel_coords, num_neighbors_Lsp = 6)
#' 
#' # Construct 26-neighbor Laplacian (corner-connected)
#' L_sparse_full <- .make_voxel_graph_laplacian_weighted(voxel_coords, num_neighbors_Lsp = 26)
#' }
#' 
#' @keywords internal
.make_voxel_graph_laplacian_weighted <- function(voxel_coords_matrix,
                                           num_neighbors_Lsp,
                                           distance_engine = c("euclidean", "ann_euclidean"),
                                           ann_threshold = 10000,
                                           weight_scheme = c("binary", "gaussian")) {
  
  # Input validation
  if (!is.matrix(voxel_coords_matrix)) {
    stop("voxel_coords_matrix must be a matrix")
  }
  
  if (ncol(voxel_coords_matrix) != 3) {
    stop("voxel_coords_matrix must have exactly 3 columns (x, y, z coordinates)")
  }
  
  if (!is.numeric(num_neighbors_Lsp) || length(num_neighbors_Lsp) != 1 || 
      num_neighbors_Lsp < 1 || num_neighbors_Lsp != round(num_neighbors_Lsp)) {
    stop("num_neighbors_Lsp must be a positive integer")
  }
  
  V <- nrow(voxel_coords_matrix)
  
  if (V < 2) {
    stop("voxel_coords_matrix must have at least 2 rows (voxels)")
  }
  
  # Ensure we don't request more neighbors than available
  k_actual <- min(num_neighbors_Lsp, V - 1)
  
  if (k_actual < num_neighbors_Lsp) {
    warning(sprintf(
      "Requested %d neighbors but only %d other voxels available. Using k=%d.",
      num_neighbors_Lsp, V - 1, k_actual
    ))
  }
  
  distance_engine <- match.arg(distance_engine)
  weight_scheme <- match.arg(weight_scheme)
  use_ann <- FALSE
  if (distance_engine == "ann_euclidean") {
    if (requireNamespace("RcppHNSW", quietly = TRUE)) {
      use_ann <- TRUE
    } else {
      warning(
        "distance_engine='ann_euclidean' requires the RcppHNSW package; falling back to exact search"
      )
    }
  } else if (distance_engine == "euclidean" && V > ann_threshold &&
             requireNamespace("RcppHNSW", quietly = TRUE)) {
    use_ann <- TRUE
  }

  # Store distances if using weighted scheme
  nn_distances <- NULL
  
  if (use_ann) {
    ann_res <- RcppHNSW::hnsw_knn(voxel_coords_matrix, k = k_actual + 1)
    nn_indices <- ann_res$idx[, -1, drop = FALSE]
    if (weight_scheme == "gaussian") {
      nn_distances <- ann_res$dist[, -1, drop = FALSE]
    }
  } else {
    # Try our C++ implementation first, fall back to RANN if needed
    if (exists("knn_search_cpp", mode = "function")) {
      nn_res <- knn_search_cpp(t(voxel_coords_matrix), t(voxel_coords_matrix), k_actual + 1)
      nn_indices <- t(nn_res$idx)[, -1, drop = FALSE]
      if (weight_scheme == "gaussian") {
        nn_distances <- t(nn_res$dist)[, -1, drop = FALSE]
      }
    } else if (requireNamespace("RANN", quietly = TRUE)) {
      nn <- RANN::nn2(voxel_coords_matrix, k = k_actual + 1)
      nn_indices <- nn$nn.idx[, -1, drop = FALSE]
      if (weight_scheme == "gaussian") {
        nn_distances <- nn$nn.dists[, -1, drop = FALSE]
      }
    } else {
      stop("No k-NN engine available: install either RcppHNSW or RANN, or ensure manifoldhrf is properly compiled.")
    }
  }
  
  # Step 2: Construct sparse adjacency matrix W
  # For undirected graph, we need to ensure symmetry
  
  # Create triplet form (i, j, value) for sparse matrix construction
  # Each voxel i is connected to its k nearest neighbors
  i_indices <- rep(seq_len(V), each = k_actual)
  j_indices <- as.vector(t(nn_indices))
  
  # Compute edge weights based on weight_scheme
  if (weight_scheme == "binary") {
    values <- rep(1, length(i_indices))
  } else if (weight_scheme == "gaussian") {
    # Use Gaussian weights: exp(-d²/σ²)
    # σ is the median distance to k-th nearest neighbor
    k_distances <- nn_distances[, k_actual]
    sigma <- median(k_distances, na.rm = TRUE)
    
    # Flatten distances
    all_distances <- as.vector(t(nn_distances))
    values <- exp(-(all_distances^2) / (sigma^2))
  }
  
  # Create initial adjacency matrix (may not be symmetric)
  W_directed <- Matrix::sparseMatrix(
    i = i_indices,
    j = j_indices,
    x = values,
    dims = c(V, V)
  )
  
  # Make symmetric by taking max(W, W^T)
  # This ensures if i is a neighbor of j OR j is a neighbor of i, they're connected
  W_directed_t <- Matrix::t(W_directed)
  W_symmetric <- pmax(W_directed, W_directed_t)
  
  # Step 3: Compute degree matrix D
  degrees <- Matrix::rowSums(W_symmetric)
  D_sparse <- Matrix::Diagonal(x = degrees)
  
  # Step 4: Compute Laplacian L = D - W
  L_sp_sparse_matrix <- D_sparse - W_symmetric
  
  # Ensure the result is stored as a sparse matrix
  L_sp_sparse_matrix <- Matrix::drop0(L_sp_sparse_matrix)  # Remove explicit zeros
  
  return(L_sp_sparse_matrix)
}

#' Apply Spatial Smoothing to Manifold Coordinates (Core)
#'
#' Spatially smooths manifold coordinates across voxels using a graph Laplacian
#' regularization approach.
#'
#' @param Xi_ident_matrix The m x V matrix of identifiability-constrained manifold 
#'   coordinates from Component 1, where m is manifold dimensionality and V is 
#'   number of voxels
#' @param L_sp_sparse_matrix The V x V sparse graph Laplacian matrix from 
#'   make_voxel_graph_laplacian_core
#' @param lambda_spatial_smooth Spatial smoothing strength parameter (scalar). 
#'   Higher values produce more smoothing. When lambda = 0, returns the 
#'   original coordinates without smoothing.
#'   
#' @return Xi_smoothed_matrix The m x V matrix of spatially smoothed manifold 
#'   coordinates
#'   
#' @details This function implements Component 2, Step 2 of the M-HRF-LSS pipeline.
#'   For each manifold dimension, it solves the regularization problem:
#'   (I + lambda * L) * xi_smooth = xi_ident
#'   This encourages nearby voxels to have similar manifold coordinates while
#'   preserving the overall structure. The identity matrix ensures the solution
#'   remains close to the original coordinates.
#'   
#' @examples
#' \dontrun{
#' # Example with synthetic data
#' m <- 5   # manifold dimensions
#' V <- 100 # voxels
#' 
#' # Create example manifold coordinates
#' Xi_ident <- matrix(rnorm(m * V), m, V)
#' 
#' # Create simple grid coordinates
#' coords <- expand.grid(x = 1:10, y = 1:10, z = 1)
#' voxel_coords <- as.matrix(coords[seq_len(V), ])
#' 
#' # Create Laplacian
#' L_sparse <- make_voxel_graph_laplacian_core(voxel_coords, num_neighbors_Lsp = 8)
#' 
#' # Apply smoothing
#' Xi_smoothed <- apply_spatial_smoothing_core(Xi_ident, L_sparse, lambda_spatial_smooth = 0.1)
#' }
#' 
#' @export
apply_spatial_smoothing_core <- function(Xi_ident_matrix,
                                        L_sp_sparse_matrix,
                                        lambda_spatial_smooth) {
  
  # Input validation
  if (!is.matrix(Xi_ident_matrix)) {
    stop("Xi_ident_matrix must be a matrix")
  }
  
  if (!inherits(L_sp_sparse_matrix, "Matrix")) {
    stop("L_sp_sparse_matrix must be a sparse matrix from the Matrix package")
  }
  
  if (!is.numeric(lambda_spatial_smooth) || length(lambda_spatial_smooth) != 1 || 
      lambda_spatial_smooth < 0) {
    stop("lambda_spatial_smooth must be a non-negative scalar")
  }
  
  # Get dimensions
  m <- nrow(Xi_ident_matrix)
  V <- ncol(Xi_ident_matrix)
  
  # Check dimension compatibility
  if (nrow(L_sp_sparse_matrix) != V || ncol(L_sp_sparse_matrix) != V) {
    stop(sprintf(
      "L_sp_sparse_matrix must be %d x %d to match Xi_ident_matrix, but is %d x %d",
      V, V, nrow(L_sp_sparse_matrix), ncol(L_sp_sparse_matrix)
    ))
  }
  
  # Early return for lambda = 0 (no smoothing)
  if (lambda_spatial_smooth == 0) {
    return(Xi_ident_matrix)
  }
  
  # Early return if there are no manifold dimensions or voxels
  if (m == 0 || V == 0) {
    return(Xi_ident_matrix)
  }
  
  # Create system matrix: (I + lambda * L)
  # I_V is the identity matrix of size V x V
  I_V <- Matrix::Diagonal(n = V)
  A_system <- I_V + lambda_spatial_smooth * L_sp_sparse_matrix
  
  # The system is symmetric positive definite, so we can use Cholesky
  # This is faster than generic solve()
  tryCatch({
    # Attempt Cholesky factorization for SPD system
    A_factor <- Matrix::Cholesky(A_system, perm = TRUE, LDL = FALSE)
    Xi_smoothed_t <- Matrix::solve(A_factor, t(Xi_ident_matrix), system = "A")
  }, error = function(e) {
    # Fall back to standard solve if Cholesky fails
    # (shouldn't happen for lambda > 0, but be safe)
    warning("Cholesky factorization failed, using standard solve: ", e$message)
    Xi_smoothed_t <- Matrix::solve(A_system, t(Xi_ident_matrix))
  })
  
  # Convert back to regular matrix and transpose to m x V
  Xi_smoothed_matrix <- t(as.matrix(Xi_smoothed_t))
  
  return(Xi_smoothed_matrix)
}
