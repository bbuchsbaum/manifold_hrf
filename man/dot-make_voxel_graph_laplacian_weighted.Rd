% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_spatial_smoothing.R
\name{.make_voxel_graph_laplacian_weighted}
\alias{.make_voxel_graph_laplacian_weighted}
\title{Construct Graph Laplacian for Voxel Coordinates with Weights (Internal)}
\usage{
.make_voxel_graph_laplacian_weighted(
  voxel_coords_matrix,
  num_neighbors_Lsp,
  distance_engine = c("euclidean", "ann_euclidean"),
  ann_threshold = 10000,
  weight_scheme = c("binary", "gaussian")
)
}
\arguments{
\item{voxel_coords_matrix}{A V x 3 matrix of voxel coordinates, where V is
the number of voxels and columns are x, y, z coordinates}

\item{num_neighbors_Lsp}{Number of nearest neighbors for graph construction
(e.g., 6 for face-connected, 18 for edge-connected, 26 for corner-connected)}

\item{distance_engine}{Choice of distance computation engine. "euclidean"
performs an exact search while "ann_euclidean" uses approximate nearest
neighbors via the \pkg{RcppHNSW} package. If \pkg{RcppHNSW} is not available,
the function falls back to the exact search and issues a warning.}

\item{ann_threshold}{When \code{distance_engine = "euclidean"}, use the
approximate search for datasets larger than this threshold if
\pkg{RcppHNSW} is installed.}

\item{weight_scheme}{Weight scheme for the adjacency matrix. "binary" uses
0/1 weights, "gaussian" uses exp(-d²/σ²) where σ is the median distance
to the k-th nearest neighbor across all voxels.}
}
\value{
L_sp_sparse_matrix A V x V sparse Laplacian matrix (Matrix::sparseMatrix)
}
\description{
Creates a sparse graph Laplacian matrix for spatial smoothing based on
k-nearest neighbors in 3D voxel space. This is an internal version that
supports weighted graphs.
}
\details{
This function implements Component 2, Step 1 of the M-HRF-LSS pipeline.
It constructs a k-nearest neighbor graph from voxel coordinates and computes
the combinatorial graph Laplacian L = D - W, where W is the symmetrized
binary adjacency matrix and D is the degree matrix. The initial k-NN graph
is directed (each voxel points to its k nearest neighbors), which is then
symmetrized by taking the maximum of W and its transpose. The resulting
Laplacian is used for spatial smoothing of manifold coordinates.
}
\examples{
\dontrun{
# Create example 3D grid of voxels
coords <- expand.grid(x = 1:10, y = 1:10, z = 1:5)
voxel_coords <- as.matrix(coords)

# Construct 6-neighbor Laplacian (face-connected)
L_sparse <- .make_voxel_graph_laplacian_weighted(voxel_coords, num_neighbors_Lsp = 6)

# Construct 26-neighbor Laplacian (corner-connected)
L_sparse_full <- .make_voxel_graph_laplacian_weighted(voxel_coords, num_neighbors_Lsp = 26)
}

}
\keyword{internal}
