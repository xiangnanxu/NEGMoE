###################################################################
# function to fitting NEGMoE with different initialization
###################################################################

.fitNEGMoE0 = function(X, Z, y, K, type, ...){

  n <- nrow(Z)
  q = ncol(Z)

  if(type == "icov"){
    res <- .EM_GMoE0(X = X, Z = Z, Y = y, K = K, ...)
    res0 <- list(Theta = res$Theta, S = res$S, mu = res$mu)
  }else if(type == "pseudoLL"){
    res <- .EM_MRMoE0(X = X, Z = Z, Y = y, K = K, ...)
    res0 <- list(W = res$W, sgm = res$sgm)
  }else if(type == "mElnet"){
    res <- .EM_MRMoE0(X = X, Z = Z, Y = y, K = K, ...)
    res0 <- list(W = res$W, sgm = res$sgm)
  }else{
    message("Please choose type from icov, pseudoLL or mElnet")
    return()
  }

  return(list(experts = res0, gamma = res$V, r_i = res$r_i, pi = res$pi,
              obj = res$obj, stop_cond = res$stop_cond, LL0 = res$LL0,
              PLL0 = res$PLL0, lambda1 = res$lambda1, lambda2 = res$lambda2,
              alpha1 = res$alpha1, alpha2 = res$alpha2))
}

###################################################################
# function to fitting NEGMoE with small EM
###################################################################
#' Fit NEGMoE model
#' @description This function fit NEGMoE model using EM algorithm.
#' @param NEGMoE_obj a NEGMoE object contain data and parameters for fitting.
#' @param gamma_init initial value of parameters in gating network.
#' @param beta_init initial values of parameters in experts network.
#' @param ... other parameters can be passed to fitNEGMoE.
#' See \code{\link{createParameterList}}
#' @return a NEGMoE object contain fitted result including parameters
#' beta a list contain coefficients of experts network in all levels.
#' gamma a matrix of coefficients in gating network.
#' r_i an array of inherit partition result in EM algorithm.
#' pen_fac estimated penalty factor for lambda1 (adapt = T). The (balanced)
#' lambda1 in NEGMoE is lambda1*pen_fac
#' latent estimated latent classes.
#'
#' @seealso \code{\link{NEGMoE_buildFromList}},
#' \code{\link{NEGMoE_buildFromPhyloseq}}, \code{\link{createParameterList}}
#' @export

fitNEGMoE <- function(NEGMoE_obj, beta_init = NULL, gamma_init = NULL,
                      ...){

  X = NEGMoE_obj@X
  Z = NEGMoE_obj@Z
  y = NEGMoE_obj@Y

  NEGMoE_result = do.call(.fitNEGMoE0,
                         c(list(X = X, Z= Z, y = y, type = NEGMoE_obj@type,
                                K = NEGMoE_obj@K, W_init = beta_init,
                                V_init = gamma_init),
                           NEGMoE_obj@params))

  NEGMoE_obj@NEGMoE_output = NEGMoE_result

  return(NEGMoE_obj)
}




