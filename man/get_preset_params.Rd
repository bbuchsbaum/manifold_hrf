% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/soundness_improvements.R
\name{get_preset_params}
\alias{get_preset_params}
\title{Preset Parameter Configurations}
\usage{
get_preset_params(
  preset = c("conservative", "balanced", "aggressive"),
  data_scale = NULL
)
}
\arguments{
\item{preset}{Character string: "conservative", "balanced", or "aggressive"}

\item{data_scale}{Optional data scale for parameter adaptation}
}
\value{
List of parameters
}
\description{
Returns preset parameter configurations for different analysis styles
}
