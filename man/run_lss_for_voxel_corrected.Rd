% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_lss.R
\name{run_lss_for_voxel_corrected}
\alias{run_lss_for_voxel_corrected}
\title{Run LSS for Single Voxel - Simple Interface}
\usage{
run_lss_for_voxel_corrected(
  y_voxel,
  X_trial_list,
  h_voxel,
  TR = 2,
  lambda = 1e-06
)
}
\arguments{
\item{y_voxel}{Data vector for single voxel (n x 1)}

\item{X_trial_list}{List of trial matrices (n x p each)}

\item{h_voxel}{HRF shape vector (p x 1)}

\item{TR}{Repetition time in seconds}

\item{lambda}{Ridge regularization parameter}
}
\value{
List with beta_trials vector
}
\description{
Simplified interface for trial-wise LSS that handles projection internally
}
