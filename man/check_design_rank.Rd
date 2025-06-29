% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostic_utils.R
\name{check_design_rank}
\alias{check_design_rank}
\title{Check Design Matrix Rank}
\usage{
check_design_rank(X_design, tol = 1e-10, remove_collinear = TRUE)
}
\arguments{
\item{X_design}{Design matrix to check}

\item{tol}{Tolerance for rank determination}

\item{remove_collinear}{Whether to remove collinear columns}
}
\value{
List with checked/cleaned design matrix and diagnostics
}
\description{
Checks rank deficiency in design matrices and removes collinear columns
}
