% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manifold_hrf_fit_class.R
\name{manifold_hrf_fit}
\alias{manifold_hrf_fit}
\title{Manifold HRF Fit Object}
\description{
S3 class representing the results of manifold HRF estimation.
This object contains the estimated HRF shapes, condition and trial amplitudes,
model fit statistics, and technical details about the manifold construction.
}
\details{
A \code{manifold_hrf_fit} object contains the following components:
\describe{
\item{amplitudes}{Condition-level amplitude estimates}
\item{trial_amplitudes}{Trial-wise amplitude estimates}
\item{hrf_shapes}{Matrix of estimated HRF shapes (timepoints x voxels)}
\item{fitted_values}{Predicted BOLD signal}
\item{residuals}{Residuals (observed - fitted)}
\item{model_specific}{List containing technical details:
\itemize{
\item manifold_coords: Manifold coordinates (Xi)
\item manifold: The manifold object
\item amplitudes_initial: Initial amplitude estimates
\item spatial_laplacian: Spatial graph Laplacian
\item convergence_info: Convergence information
}
}
\item{call}{The matched call}
\item{control}{Control parameters used}
\item{data_info}{Data dimensions and metadata}
\item{qc_metrics}{Optional quality control metrics}
}
}
\seealso{
\code{\link{estimate_hrf_manifold}}, \code{\link{manifold_control}}
}
