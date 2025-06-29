% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostic_utils.R
\name{monitor_condition_numbers}
\alias{monitor_condition_numbers}
\title{Monitor Condition Numbers}
\usage{
monitor_condition_numbers(
  matrix_list,
  warn_threshold = 1e+08,
  error_threshold = 1e+12
)
}
\arguments{
\item{matrix_list}{Named list of matrices to check}

\item{warn_threshold}{Threshold for warning (default 1e8)}

\item{error_threshold}{Threshold for error (default 1e12)}
}
\value{
Data frame with condition monitoring results
}
\description{
Monitors condition numbers throughout the pipeline
}
