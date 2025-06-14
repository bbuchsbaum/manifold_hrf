% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostic_utils.R
\name{track_convergence_metrics}
\alias{track_convergence_metrics}
\title{Track Convergence Metrics}
\usage{
track_convergence_metrics(
  current_values,
  previous_values = NULL,
  iteration = 1,
  metric_name = "parameters",
  history = NULL
)
}
\arguments{
\item{current_values}{Current parameter values (vector or matrix)}

\item{previous_values}{Previous parameter values}

\item{iteration}{Current iteration number}

\item{metric_name}{Name of the metric being tracked}

\item{history}{Optional existing history to append to}
}
\value{
Updated convergence history
}
\description{
Tracks convergence metrics across iterations for diagnostics
}
