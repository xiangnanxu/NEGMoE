% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GLasso.R
\name{soft_thresh}
\alias{soft_thresh}
\title{soft threshold function}
\usage{
soft_thresh(x, lambda)
}
\arguments{
\item{x}{input of vector for soft-threshold function}

\item{lambda}{penalty parameters for soft-threshold function.}
}
\value{
output is sign(x)*max(x-lambda,0)
}
\description{
soft threshold function
}
