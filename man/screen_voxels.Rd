% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robust_spatial_outlier.R
\name{screen_voxels}
\alias{screen_voxels}
\title{Screen Voxels for Quality}
\usage{
screen_voxels(Y_data, min_variance = 1e-06, max_spike_fraction = 0.1)
}
\arguments{
\item{Y_data}{n x V data matrix}

\item{min_variance}{Minimum temporal variance threshold}

\item{max_spike_fraction}{Maximum fraction of spike-like values}
}
\value{
List with keep/flag indicators and quality metrics
}
\description{
Identifies voxels that should be excluded or flagged
}
