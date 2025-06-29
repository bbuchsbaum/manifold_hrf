% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config_helpers.R
\name{suggest_parameters}
\alias{suggest_parameters}
\title{Adaptive Parameter Selection}
\usage{
suggest_parameters(Y_data, X_design_list = NULL, voxel_coords = NULL)
}
\arguments{
\item{Y_data}{n x V data matrix}

\item{X_design_list}{List of design matrices}

\item{voxel_coords}{V x 3 coordinate matrix (optional)}
}
\value{
List of suggested parameters
}
\description{
Suggests parameters based on data characteristics
}
