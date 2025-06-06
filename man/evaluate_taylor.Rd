% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taylor_utils.R
\name{evaluate_taylor}
\alias{evaluate_taylor}
\title{Evaluate a Taylor series at given points}
\usage{
evaluate_taylor(coeff, x0, x)
}
\arguments{
\item{coeff}{Numeric vector of Taylor coefficients.}

\item{x0}{Expansion point used for the coefficients.}

\item{x}{Numeric vector of evaluation points.}
}
\value{
Numeric vector of the Taylor polynomial evaluated at \code{x}.
}
\description{
Evaluate a Taylor series at given points
}
