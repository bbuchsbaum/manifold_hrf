% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_voxelfit_engine.R, R/core_voxelwise_fit.R
\name{extract_xi_beta_raw_svd_core}
\alias{extract_xi_beta_raw_svd_core}
\title{Extract raw manifold coordinates and condition amplitudes via SVD}
\usage{
extract_xi_beta_raw_svd_core(
  Gamma_coeffs_matrix,
  m_manifold_dim,
  k_conditions,
  n_jobs = 1,
  log_fun = NULL
)

extract_xi_beta_raw_svd_core(
  Gamma_coeffs_matrix,
  m_manifold_dim,
  k_conditions,
  n_jobs = 1,
  log_fun = NULL
)
}
\arguments{
\item{Gamma_coeffs_matrix}{The (km) x V coefficient matrix from
solve_glm_for_gamma_core, where rows are organized by condition and
manifold dimension}

\item{m_manifold_dim}{Manifold dimensionality (m)}

\item{k_conditions}{Number of conditions (k)}

\item{n_jobs}{Number of parallel jobs for voxel processing (default 1).}

\item{log_fun}{Optional logging function to report quality metrics. If NULL,
no logging is performed.}
}
\value{
list with Xi_raw_matrix (m x V) and Beta_raw_matrix (k x V)

A list containing:
\itemize{
\item \code{Xi_raw_matrix}: m x V matrix of raw manifold coordinates
\item \code{Beta_raw_matrix}: k x V matrix of raw condition amplitudes
\item \code{quality_metrics}: Diagnostics from the robust SVD
}
}
\description{
Extracts raw manifold coordinates (Xi) and condition amplitudes (Beta) from
gamma coefficients using singular value decomposition.
}
\details{
This function implements Component 1, Step 4 of the M-HRF-LSS pipeline.
It now delegates the SVD step to \code{extract_xi_beta_raw_svd_robust}, which
adds condition number checks and fallback strategies. The first singular
value and associated vectors are used to decompose gamma into Xi (HRF shape)
and Beta (amplitude) components. Near-zero singular values are handled by
setting the corresponding Xi and Beta values to zero. Quality metrics from
the robust SVD can optionally be logged via \code{log_fun}.
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
