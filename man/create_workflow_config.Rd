% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config_helpers.R
\name{create_workflow_config}
\alias{create_workflow_config}
\title{Create Full Workflow Configuration}
\usage{
create_workflow_config(
  preset = "balanced",
  custom_params = NULL,
  data_checks = TRUE,
  output_dir = NULL
)
}
\arguments{
\item{preset}{Base preset to use}

\item{custom_params}{List of custom parameter overrides}

\item{data_checks}{Whether to enable data validation checks}

\item{output_dir}{Directory for outputs (QC reports, logs, etc.)}
}
\value{
Complete workflow configuration
}
\description{
Creates a complete workflow configuration with all settings
}
