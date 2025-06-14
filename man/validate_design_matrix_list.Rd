% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_lss.R, R/input_validation.R
\name{validate_design_matrix_list}
\alias{validate_design_matrix_list}
\title{Validate design matrix list}
\usage{
validate_design_matrix_list(X_list, n_timepoints = NULL)

validate_design_matrix_list(X_list, n_timepoints = NULL)
}
\arguments{
\item{X_list}{List of design matrices}

\item{n_timepoints}{Expected number of rows for each matrix (optional)}
}
\value{
List with dimension information
}
\description{
Ensures that each design matrix is numeric with consistent dimensions and
aligned with the expected number of timepoints.
}
\keyword{internal}
