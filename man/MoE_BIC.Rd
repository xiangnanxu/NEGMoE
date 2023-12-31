% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Criterion.R
\name{MoE_BIC}
\alias{MoE_BIC}
\title{Calculate BIC of a RMoE result}
\usage{
MoE_BIC(EM_res, eps = 1e-05)
}
\arguments{
\item{EM_res}{A list of RMoE result.}

\item{eps}{A threshold for minimal coefficient for calculate degree of freedom.}
}
\value{
A numerical variable of BIC statistics of the fitted RMoE model.
}
\description{
Calculate BIC of a RMoE result
}
