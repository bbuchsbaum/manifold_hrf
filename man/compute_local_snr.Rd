% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robust_spatial_outlier.R
\name{compute_local_snr}
\alias{compute_local_snr}
\title{Compute Local SNR for Adaptive Smoothing}
\usage{
compute_local_snr(
  Y_data,
  Y_predicted = NULL,
  method = c("temporal_variance", "residual")
)
}
\arguments{
\item{Y_data}{n x V data matrix}

\item{Y_predicted}{n x V predicted data (optional)}

\item{method}{Method for SNR estimation}
}
\value{
Vector of SNR values for each voxel
}
\description{
Estimates signal-to-noise ratio for each voxel to guide smoothing strength
}
\keyword{internal}
