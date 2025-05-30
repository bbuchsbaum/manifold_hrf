% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fallback_memory_safety.R
\name{check_and_scale_data}
\alias{check_and_scale_data}
\title{Check Data Scaling and Adjust if Needed}
\usage{
check_and_scale_data(Y_data, target_scale = 100)
}
\arguments{
\item{Y_data}{Data matrix}

\item{target_scale}{Target scale for data}
}
\value{
List with scaled data and scaling factor
}
\description{
Ensures data is in reasonable range for numerical stability
}
