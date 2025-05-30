% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robust_voxelwise_fit.R
\name{smart_initialize}
\alias{smart_initialize}
\title{Smart Initialization for Voxelwise Fit}
\usage{
smart_initialize(
  Y_data,
  X_condition_list,
  hrf_canonical,
  use_spatial_clusters = TRUE,
  voxel_coords = NULL,
  m_manifold_dim = 5
)
}
\arguments{
\item{Y_data}{n x V data matrix}

\item{X_condition_list}{List of condition design matrices}

\item{hrf_canonical}{Canonical HRF for initialization}

\item{use_spatial_clusters}{Whether to use spatial clustering}

\item{voxel_coords}{Voxel coordinates for clustering}

\item{m_manifold_dim}{Manifold dimension for Xi initialization}
}
\value{
List with initial Xi and Beta matrices
}
\description{
Provides intelligent starting values for the manifold fit
}
