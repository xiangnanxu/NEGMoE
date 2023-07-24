#' check validity
#' @param object A NEGMoE object.
check_validity <- function(object){
  
  n_M = nrow(object@X)
  n_N = nrow(object@Z)
  n_R = length(object@Y)
  
  if(!all(n_M == n_N)){
    errors <- c("Numbers of samples in X not equal
                to Z!")
    return(errors)
  }
  
  if(!is.null(object@taxLevel)){
    if(length(object@X) != length(object@taxLevel)){
      errors <- c("Selected taxa level should be the same as
                  names in microbiome list!")
      return(errors)
    }
  }
  if(is.integer(object@K)){
    errors = c("Number of latent classes should be integer!")
  }
  
  return(TRUE)
}

#' @importClassesFrom S4Vectors DataFrame DFrame
setClassUnion("data.frameORNULL", c("data.frame", "NULL"))

setClassUnion("characterORNULL", c("character", "NULL"))

setClassUnion("listORNULL", c("list", "NULL"))

setClassUnion("matrixORNULL", c("matrix", "data.frame", "vector", "NULL"))

#' NEGMoE class
#' @slot X A matrix of microbiome data, rows are samples.
#' @slot Z A matrix of nutrition intake.
#' @slot Y A vector of response or NULL.
#' @slot type Type of model, can be chosen from "icov", "pseudoLL" or "mElnet".
#' "mElnet" corresponding to multivariate elastic net model. (Y must given)
#' "icov" corresponding to sparse inverse covariance matrix estimation. (graphical lasso)
#' "pseudoLL" corresponding to pseudo likelihood method. (mb in huge package)
#' By default, will use "pseudoLL".
#' @slot params A list of parameters for fitting NEGMoE.
#' @slot K A number of fitting componenets in NEGMoE (>=2).
#' @slot cvParams A list of parameters for selecting parameters of NEGMoE
#'  using cross validation.
#' @slot cvResult A dataframe of cross validation results of different parameters
#'  of fitting NEGMoE.
#' @slot NEGMoE_output A list NEGMoE fitting result.
#' @slot taxLevel A character indicates selected taxa level.
#' @slot ResponseLevel A character indicates levels of response variables.
#' @slot standardize A Logical variable indicate whether standardize input data.
#' @importFrom methods new
NEGMoE <- setClass("NEGMoE",
                  slots = c(X = "matrix",
                            Z = "matrix",
                            Y = "matrixORNULL",
                            type = "character",
                            params = "listORNULL",
                            K = "numeric",
                            cvParams = "listORNULL",
                            cvResult = "listORNULL",
                            NEGMoE_output = "list",
                            taxLevel = "characterORNULL",
                            taxTab = "data.frameORNULL",
                            standardize = "logical",
                            .transformation = "listORNULL"
                  ),
                  prototype = list(
                    params = list(),
                    cvParams = list(),
                    .transformation = list(method = "none")
                  ))

#' @importFrom S4Vectors coolcat
#'

setMethod("show", "NEGMoE", function(object) {
  
#  L = length(object@X)
  X_char = paste0("X contains number of variables: ", dim(object@X)[2])

  X_char = paste0(X_char, "\n")
  
  cat("Numer of samples: ", nrow(object@X), '\n')
  cat(X_char)
  cat("Z dim:", dim(object@Z), "\n")
  cat("Number of latent class K: ", object@K, '\n')
  S4Vectors::coolcat("NEGMoE output names (%d): %s\n",
                     names(object@NEGMoE_output))
  
  S4Vectors::coolcat("param names (%d): %s\n", names(object@params))
  S4Vectors::coolcat("cross validation param names (%d): %s\n",
                     names(object@cvParams))
})

#' @importFrom S4Vectors setValidity2
#'
S4Vectors::setValidity2("NEGMoE", check_validity)

