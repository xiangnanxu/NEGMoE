% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitNEGMoE.R
\name{fitNEGMoE}
\alias{fitNEGMoE}
\title{Fit NEGMoE model}
\usage{
fitNEGMoE(NEGMoE_obj, beta_init = NULL, gamma_init = NULL, ...)
}
\arguments{
\item{NEGMoE_obj}{a NEGMoE object contain data and parameters for fitting.}

\item{beta_init}{initial values of parameters in experts network.}

\item{gamma_init}{initial value of parameters in gating network.}

\item{...}{other parameters can be passed to fitNEGMoE.
See \code{\link{createParameterList}}}
}
\value{
a NEGMoE object contain fitted result including parameters
beta a list contain coefficients of experts network in all levels.
gamma a matrix of coefficients in gating network.
r_i an array of inherit partition result in EM algorithm.
pen_fac estimated penalty factor for lambda1 (adapt = T). The (balanced)
lambda1 in NEGMoE is lambda1*pen_fac
latent estimated latent classes.
}
\description{
This function fit NEGMoE model using EM algorithm.
}
\seealso{
\code{\link{NEGMoE_buildFromList}},
\code{\link{NEGMoE_buildFromPhyloseq}}, \code{\link{createParameterList}}
}
