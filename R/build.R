###################################################################
# function to build NEGMoE object from List, dataframe or matrix
###################################################################
#' Build NEGMoE object from list, dataframe or matrix
#' @description This function build NEGMoE object from List, dataframe or matrix.
#' @param X A list of input matrix (Microbiome)
#' data of different taxonomic level. If it is a dataframe or matrix, will
#' convert to a list automatically.
#' @param Z A dataframe of matrix of gating networks.
#' @param Y A multivariate outcome matrix. If is NULL, type will be automatically change to "pseudoLL".
#' @param type Type of model, can be chosen from "icov", "pseudoLL" or "mElnet".
#' "mElnet" corresponding to multivariate elastic net model. (Y must given)
#' "icov" corresponding to sparse inverse covariance matrix estimation. (graphical lasso)
#' "pseudoLL" corresponding to pseudo likelihood method. (mb in huge package)
#' By default, will use "pseudoLL".
#' @param K A integer of number of latent class.
#' @param cvParams A list of cross validation parameters.
#' @param taxLevel A character of selected name of taxonomic levels.
#' @param taxTab A dataframe of taxonomic table.
#' @param standardize Logical flag for x variable standardization.
#' Default is standardize=TRUE.
#' @param ... Other parameters can pass to NEMoE_buildFromList.
#' See \code{\link{createParameterList}}
#' @return A NEGMoE object.
#' @export
#'
NEGMoE_buildFromList <- function(X, Z, type = "pseudoLL", Y = NULL,
                                K = NULL, cvParams = list(),
                                taxLevel = NULL, taxTab = NULL,
                                standardize = TRUE, ...){

  if(is.data.frame(X) || is.matrix(X)){
    X = as.matrix(X)
  }

  if(is.data.frame(Z) || is.matrix(Z)){
    Z = as.matrix(Z)
  }

  if(!is.null(Y)){
    if(is.vector(Y)){
      Y <- as.matrix(as.numeric(Y))
    }else{
      Y <- as.matrix(Y)
    }
    message("Y is given, will use mElnet mode")
    type = "mElnet"
  }

  if(!(type %in% c("icov","pseudoLL"))){
    message("type must be chosen from icov, pseudoLL, or mElnet.
            Automatically change to pseudoLL")
    type = "pseudoLL"
  }else if((type == "mElnet")&(is.null(Y))){
    message("type is mElnet but Y not provided. Automatically change to pseudoLL")
    type = "pseudoLL"
  }

  params = createParameterList(...)

  if(!length(cvParams)){
    cvParams = createCVList()
  }

  if(is.null(K)){
    K = 2
  }

#  id0 <- list()
#  for(i in 1:L){
#    M_temp <- filterComp(Microbiome[[i]], var, 1e-9, TRUE)
#    Microbiome[[i]] <- M_temp$X
#    id0[[i]] = as.logical(M_temp$id)
#    names(id0[[i]]) <- colnames(Microbiome[[i]])
#  }

  if(standardize){

    .transformation = list(method = "scale", mu_X = c(), sd_X = c(),
                           mu_Z = c(), sd_Z = c())

    Z <- scale(Z)
    .transformation$mu_Z <- attr(Z,"scaled:center")
    .transformation$sd_Z <- attr(Z,"scaled:scale")

#    for(i in 1:L){
      X <- scale(X)
      .transformation$mu_X = attr(X, "scaled:center")
      .transformation$sd_X = attr(X, "scaled:scale")
#    }
#    .transformation$keepid = id0

  }else{
    .transformation = list(method = "none", mu_X = c(), sd_X = c(),
                           mu_Z = c(), sd_Z = c())
  }

  NEGMoE_obj = NEGMoE(X = X, Z = Z, Y = Y, type = type, params = params,
                      cvParams = cvParams, K = K, taxLevel = taxLevel,
                      taxTab = taxTab, .transformation = .transformation,
                      standardize = standardize)

  return(NEGMoE_obj)

}

###################################################################
# function to create list of parameters for fitting NEMoE
###################################################################
#' Create list of parameters for fitting NEMoE
#' @description This function create parameters that put into fitNEMoE function.
#' @param lambda1 Penalty regularizer for the experts.
#' @param lambda2 Penalty regularizer for the gating network.
#' @param alpha1 Elastic net penalty value for experts.
#' @param alpha2 Elastic net penalty value for gating network.
#' @param stop_all Method of stop criterion. If stop_all = TRUE means that either
#'  coefficient or loss function converge. If stop_all = FALSE means that both
#'  coefficient and loss function converge.
#' @param EM_alg Method for Expecation maximization update.
#'  Can be chosen from "EM", "CEM", "GEM", "SEM", "SAEM", "GEM".
#'  By default is "EM".
#' @param itmax Maximium number of iteration in fitting NEMoE. By default is 100.
#' @param itmin Minimium number of iteration in fitting NEMoE. By default is 3.
#' @param backtracking Whether use backtracking during the iteration.
#' Can be chosen from
#' @param init Method for initialization.
#' Can be chosen from "rand", "kmeans" and "glmnet".
#' If init="rand" will use a dirichlet distribution initialing the latent class.
#' If init="kmeans" the latent class will initialized using kmeans
#' clustering of input for gating network.
#' If init="glmnet" will use a equal probability
#' with lasso as its corresponding coefficients in experts network.
#' @param beta_max Maximal of coefficient. By default is 10.
#' @param adapt whether use adaptive mode in optimization. By default is TRUE.
#' @param verbose A logical input indicating whether the intermediate
#' steps will be printed.
#' @param early_stop A logical input indicate whether early stop when one of
#' the fitted latent class have zero variables selected (to save time).
#' By default is TRUE.
#' @return A list contain parameters in fitting NEMoE.
#' @export
#' @seealso \code{\link{NEMoE_buildFromList}},
#' \code{\link{NEMoE_buildFromPhyloseq}}, \code{\link{fitNEMoE}}

createParameterList <- function(lambda1 = 0.2, lambda2 = 0.03, alpha1 = 1,
                                alpha2 = 1, stop = "all", gamma_func = inv_func, beta_max = 10,
                                solver = "glmnet", stop_eps1 = 1e-2, stop_eps2 = 1e-2,
                                itmax = 1e2, itmin = 3, itmax0 = 10, itmax1 = 1e2, itmax2 = 1e2,
                                EM_opt = "EM", backtracking = "comp",adapt = T,
                                init = "kmeans",ratio_max = 20, method = "icov",
                                thresh_cor = 0.3, thresh_eig = 0.1, verbose = TRUE){

  params = list(lambda1 = lambda1, lambda2 = lambda2,
                alpha1 = alpha1, alpha2 = alpha2,
                stop = stop, gamma_func = inv_func, beta_max = beta_max,
                solver = solver, stop_eps1 = stop_eps1, stop_eps2 = stop_eps2,
                itmax = itmax, itmin = itmin, itmax0 = itmax0, itmax1 = itmax1,
                itmax2 = itmax2, EM_opt = EM_opt, backtracking = backtracking,
                adapt = adapt, init = init, ratio_max = ratio_max, method = method,
                thresh_cor = thresh_cor, thresh_eig = thresh_eig, verbose = verbose)

  return(params)
}

###################################################################
# function to create list of parameters for cross validate NEMoE
###################################################################
#' Create list of parameters for cross validate NEMoE
#' @description This function create parameters that put into cvNEMoE function.
#' @param g1 Numbers of parameters lambda1 for validation.
#' @param lambda2_seq A vector of candidates of lambda2.
#' @param shrink A number of shrinkage of selected lambda1. By default is 0.5.
#' @param crit_eval A vector of method for Evaluation metric of NEMoE object.
#'  Can be chosen from statistics "AIC", "BIC", "ICL1", "ICL2", "eBIC",
#'  "mAIC", "mBIC", "mICL1", "mICL2" and cross validation result "accuracy",
#'  "D.square", "TPR", "TNR", "F1" and "auc".
#'  If method = "all", all of the evaluation metric will be use.
#'  By default is "all.
#' @param crit_sel Method use for select parameters.
#' Should be a member in crit_eval.
#' @param track Whether output the result of choosing parameters lambda1.
#' @param itmax_lambda Maximal iterations of generation lambda1.
#' By default is 3.
#' @param itmax_cv Maximal iterations of calculate cross validation metric.
#' By default is 20.
#' @param itmax_fit Maximal iterations of fitting NEMoE inside the evaluation.
#' By default is 50.
#' @return A list contain parameters in fitting NEMoE.
#' @export
#' @seealso \code{\link{NEMoE_buildFromList}},
#' \code{\link{NEMoE_buildFromPhyloseq}}, \code{\link{fitNEMoE}},
#' \code{\link{calcCriterion}}, \code{\link{cvNEMoE}}

createCVList <- function(g1 = 10,
                         lambda2_seq = c(0.005, 0.014, 0.016, 0.023, 0.025),
                         shrink = 0.5, crit_eval = "all", crit_sel = "outLL",
                         track = FALSE, itmax_lambda = 3,
                         itmax_cv = 20, itmax_fit = 50, prob_thresh = 0.5){

  cvParams = list(g1 = g1, lambda2_seq = lambda2_seq, shrink = shrink,
                  crit_eval = crit_eval, crit_sel = crit_sel,
                  track = track, itmax_lambda = itmax_lambda,
                  itmax_cv = itmax_cv, itmax_fit = itmax_fit)
  return(cvParams)
}



###################################################################
# function to set parameters list of NEMoE object
###################################################################
#' Set parameters of fitting NEMoE object
#' @description This function set parameters in NEMoE object.
#' @param NEMoE_obj A NEMoE object.
#' @param ... Other parameters can pass to setParam.
#'
#' @return A NEMoE with user input parameters.
#' @export
#' @seealso \code{\link{createParameterList}}

setParam <- function(NEGMoE_obj, ...){

  params <- createParameterList(...)
  NEGMoE_obj@params <- params

  return(NEGMoE_obj)

}

###################################################################
# function to get fitted result of NEGMoE object
###################################################################
#' Get coefficients of fitted NEGMoE object
#' @description This function get coefficients in NEGMoE object.
#' @param NEMoE_obj A NEGMoE object.
#'
#' @return A list contain fitted result.
#' @export
#' @seealso \code{\link{fitNEMoE}}

getCoef <- function(NEGMoE_obj){

  gating <- NEGMoE_obj@NEMoE_output$gamma
  experts <- NEGMoE_obj@NEMoE_output$experts

  return(list(coef.gating = gating, coef.experts = experts))

}


###################################################################
# function to get fitted Likelihood of NEMoE object
###################################################################
#' Get Likelihood of fitted NEMoE object
#' @description This function get Likelihood result in NEMoE object.
#' @param NEMoE_obj A NEMoE object.
#'
#' @return A dataframe contain fitted result.
#' @export
#' @seealso \code{\link{fitNEMoE}}

getLL <- function(NEGMoE_obj){

  return(NEGMoE_obj@NEGMoE_output$obj)

}
