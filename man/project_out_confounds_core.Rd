% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_voxelfit_engine.R, R/core_voxelwise_fit.R
\name{project_out_confounds_core}
\alias{project_out_confounds_core}
\title{Project out confound regressors from data and design matrices}
\usage{
project_out_confounds_core(
  Y_data_matrix,
  X_list_of_matrices,
  Z_confounds_matrix = NULL
)

project_out_confounds_core(
  Y_data_matrix,
  X_list_of_matrices,
  Z_confounds_matrix = NULL
)
}
\arguments{
\item{Y_data_matrix}{An n x V matrix of BOLD data (n timepoints, V voxels)}

\item{X_list_of_matrices}{A list of k matrices, each n x p (design matrices
for k conditions)}

\item{Z_confounds_matrix}{An n x q_confound matrix of confound regressors,
or NULL if no confounds}
}
\value{
list with Y_proj_matrix and X_list_proj_matrices

A list containing:
\itemize{
\item \code{Y_proj_matrix}: The n x V projected data matrix
\item \code{X_list_proj_matrices}: List of k projected design matrices
}
}
\description{
Projects out confound variables from both the data matrix and design matrices
using QR decomposition.
}
\details{
This function implements Component 1, Step 1 of the M-HRF-LSS pipeline.
It uses QR decomposition to efficiently project out confound regressors from
both the data and design matrices. If \code{Z_confounds_matrix} is \code{NULL},
the original matrices are returned unchanged. Rank-deficient confound
matrices are automatically reduced to their independent columns with a
warning. Missing values are not allowed.
}
\examples{
\dontrun{
# Create example data
n <- 200  # timepoints
V <- 100  # voxels
p <- 30   # HRF length
k <- 3    # conditions

Y_data <- matrix(rnorm(n * V), n, V)
X_list <- lapply(1:k, function(i) matrix(rnorm(n * p), n, p))
Z_confounds <- cbind(1, poly(1:n, degree = 3))  # intercept + polynomial trends

# Project out confounds
result <- project_out_confounds_core(Y_data, X_list, Z_confounds)
}

}
