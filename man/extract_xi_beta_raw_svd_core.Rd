% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_voxelwise_fit.R
\name{extract_xi_beta_raw_svd_core}
\alias{extract_xi_beta_raw_svd_core}
\title{Extract Xi and Beta via SVD (Core)}
\usage{
extract_xi_beta_raw_svd_core(Gamma_coeffs_matrix, m_manifold_dim, k_conditions)
}
\arguments{
\item{Gamma_coeffs_matrix}{The (km) x V coefficient matrix from
solve_glm_for_gamma_core, where rows are organized by condition and
manifold dimension}

\item{m_manifold_dim}{Manifold dimensionality (m)}

\item{k_conditions}{Number of conditions (k)}
}
\value{
A list containing:
\itemize{
\item \code{Xi_raw_matrix}: m x V matrix of raw manifold coordinates
\item \code{Beta_raw_matrix}: k x V matrix of raw condition amplitudes
}
}
\description{
Extracts raw manifold coordinates (Xi) and condition amplitudes (Beta) from
gamma coefficients using singular value decomposition.
}
\details{
This function implements Component 1, Step 4 of the M-HRF-LSS pipeline.
For each voxel, it reshapes the gamma coefficients into an m x k matrix and
performs SVD to extract the dominant pattern. The first singular value and
vectors are used to decompose gamma into Xi (HRF shape) and Beta (amplitude)
components. Near-zero singular values are handled by setting the corresponding
Xi and Beta values to zero.
}
\examples{
\dontrun{
# Create example gamma coefficients
m <- 5   # manifold dimensions
k <- 3   # conditions
V <- 100 # voxels

# Gamma from GLM solve
gamma <- matrix(rnorm((k * m) * V), k * m, V)

# Extract Xi and Beta
result <- extract_xi_beta_raw_svd_core(gamma, m, k)
# result$Xi_raw_matrix is m x V
# result$Beta_raw_matrix is k x V
}

}
