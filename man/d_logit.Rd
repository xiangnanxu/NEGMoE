% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{d_logit}
\alias{d_logit}
\title{Derivative of logit function}
\usage{
d_logit(x, beta)
}
\arguments{
\item{x}{A matrix of input for logistic regression}

\item{beta}{A vector of coefficients.}
}
\value{
\frac{d(logit(x\beta))}{d\beta}
}
\description{
Derivative of logit function
}
