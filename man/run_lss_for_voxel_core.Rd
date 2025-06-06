% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_lss_deprecated.R
\name{run_lss_for_voxel_core}
\alias{run_lss_for_voxel_core}
\title{Run LSS for Single Voxel (Core)
Run LSS for Single Voxel (Core)}
\usage{
run_lss_for_voxel_core(
  Y_proj_voxel_vector,
  X_trial_onset_list_of_matrices,
  H_shape_voxel_vector,
  A_lss_fixed_matrix,
  P_lss_matrix,
  p_lss_vector
)
}
\arguments{
\item{Y_proj_voxel_vector}{An n x 1 vector of projected BOLD data for one voxel}

\item{X_trial_onset_list_of_matrices}{A list of T matrices, where each is an
n x p Toeplitz design matrix for a single trial}

\item{H_shape_voxel_vector}{A p x 1 vector representing the HRF shape for
this voxel}

\item{A_lss_fixed_matrix}{The n x q_lss matrix of fixed regressors}

\item{P_lss_matrix}{The q_lss x n precomputed projection matrix from
prepare_lss_fixed_components_core}

\item{p_lss_vector}{The n x 1 precomputed intercept projection vector from
prepare_lss_fixed_components_core}
}
\value{
beta_trial_voxel_vector A T x 1 vector of trial-wise beta estimates
}
\description{
Performs Least Squares Separate (LSS) estimation for a single voxel using
the Woodbury matrix identity for computational efficiency.
}
\details{
This function implements Component 3, Step 3 of the M-HRF-LSS pipeline.
For each trial, it uses the Woodbury matrix identity to efficiently compute
the beta estimate when that trial's regressor is added to the fixed model.
This avoids repeated matrix inversions and makes single-trial estimation
computationally feasible.

The Woodbury identity allows us to update the inverse when adding a single
column c to the design matrix:
(A'A + cc')^(-1) = (A'A)^(-1) - (A'A)^(-1)cc'(A'A)^(-1) / (1 + c'(A'A)^(-1)c)
}
\examples{
\dontrun{
# Setup for single voxel
n <- 200  # timepoints
p <- 30   # HRF length
T <- 50   # trials

# Data for one voxel
Y_voxel <- rnorm(n)

# Trial onsets (simplified)
X_trials <- lapply(1:T, function(t) {
  X <- matrix(0, n, p)
  # Put onset at different times
  onset_time <- 10 + (t-1) * 3
  if (onset_time + p <= n) {
    X[onset_time:(onset_time+p-1), ] <- diag(p)
  }
  X
})

# HRF shape for this voxel
H_voxel <- dgamma(0:(p-1), shape = 6, rate = 1)
H_voxel <- H_voxel / max(H_voxel)

# Fixed regressors
A_fixed <- cbind(1, seq_len(n)/n)
lss_prep <- prepare_lss_fixed_components_core(A_fixed, 1, 1e-6)

# Run LSS
betas <- run_lss_for_voxel_core(
  Y_voxel, X_trials, H_voxel, A_fixed,
  lss_prep$P_lss_matrix, lss_prep$p_lss_vector
)
}

}
