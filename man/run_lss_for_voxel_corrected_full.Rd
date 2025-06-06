% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_lss.R
\name{run_lss_for_voxel_corrected_full}
\alias{run_lss_for_voxel_corrected_full}
\title{Run LSS for Single Voxel - Full Interface}
\usage{
run_lss_for_voxel_corrected_full(
  Y_proj_voxel_vector,
  X_trial_onset_list_of_matrices,
  H_shape_voxel_vector,
  A_lss_fixed_matrix,
  P_lss_matrix,
  p_lss_vector
)
}
\arguments{
\item{Y_proj_voxel_vector}{Projected data (n x 1) with confounds removed}

\item{X_trial_onset_list_of_matrices}{List of unprojected trial matrices (n x p each)}

\item{H_shape_voxel_vector}{HRF shape (p x 1)}

\item{A_lss_fixed_matrix}{Matrix of fixed regressors used during projection (n x q)}

\item{P_lss_matrix}{Precomputed matrix from \code{prepare_lss_fixed_components_core} (q x n)}

\item{p_lss_vector}{Precomputed intercept projection vector of length n}
}
\value{
Vector of trial-wise betas
}
\description{
Complete interface using precomputed LSS components
}
