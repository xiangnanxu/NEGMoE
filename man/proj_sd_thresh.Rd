% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GLasso.R
\name{proj_sd_thresh}
\alias{proj_sd_thresh}
\title{Project a matrix to correlation space with non-zero correlation larger than thresh_cor.}
\usage{
proj_sd_thresh(X, thresh_cor = 0.3, thresh_eig = 0.1)
}
\arguments{
\item{X}{Input of the matrix}

\item{thresh_cor}{threshold for minimal non-zero correlation of output semi-definite matrix.}

\item{thresh_eig}{threshold for eigenvalue of output semi-definite matrix.}
}
\value{
A semi-definite matrix by iterative projection.
}
\description{
Project a matrix to correlation space with non-zero correlation larger than thresh_cor.
}
