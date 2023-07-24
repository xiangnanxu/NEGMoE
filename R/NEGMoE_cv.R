.lambda1Generator <- function(NEGMoE_obj, lambda2, itmax = 5, lambda1_max = 3){

  p <- ncol(NEGMoE_obj@X)
  q <- ncol(NEGMoE_obj@Z)
  y <- NEGMoE_obj@Y
  K <- NEGMoE_obj@K
  NEGMoE_temp = NEGMoE_obj
  NEGMoE_temp@params$lambda2 = lambda2
  gamma_init = NEGMoE_obj@NEGMoE_output$gamma
  NEGMoE_temp@params$verbose = F

  type = NEGMoE_obj@type
  lambda1_max = lambda1_max
  repeat{
    NEGMoE_temp@X = NEGMoE_obj@X
    NEGMoE_temp@params$lambda1 = lambda1_max
    NEGMoE_temp@params$itmax = itmax
    if(type == "icov"){
      Theta_0 = array(0, dim = c(p,p,K))
      for(i in 1:K){
        Theta_0[,,i] = diag(rep(1,p))
      }
      beta_init = Theta_0
    }else{
      beta_init <- array(0, dim = c(p + 1,p,K))
    }
    NEGMoE_temp = fitNEGMoE(NEGMoE_temp, beta_init = beta_init,
                            gamma_init = gamma_init)

    df_beta = calcdf(NEGMoE_temp)$df_beta

    if(!df_beta){
      lambda1_max = lambda1_max*0.8
    }else{
      break()
    }
  }
  return(lambda1_max)

}

###################################################################
# function to calculate evaluation metric of NEGMoE object
###################################################################
#' Evaluate fitted NEGMoE obejct
#' @description This function calculate evaluation of NEGMoE object.
#' @param NEGMoE a NEGMoE of object with fitted parameters.
#' @param crit_list a list of evaluation metric.
#' @param ... other parameters can be put into calcCriterion.
#' @return A vector of evaluation result.
#' @export
#' @seealso \code{\link{createCVList}}


calcCriterion <- function(NEGMoE_obj, crit_list= c("AIC", "BIC"),
                          ...){

  f_AIC <- .NEGMoE_AIC
  f_BIC <- .NEGMoE_BIC
  f_eBIC <- purrr::partial(.NEGMoE_BIC, type = "eBIC")
  f_mAIC <- purrr::partial(.NEGMoE_AIC, modified = TRUE)
  f_mBIC <- purrr::partial(.NEGMoE_BIC, modified = TRUE)


  if("all" %in% crit_list){
    crit_list <- c("AIC", "BIC", "eBIC", "mAIC", "mBIC", "outLL")
  }

  f_list = list(f_AIC, f_BIC, f_eBIC, f_mAIC, f_mBIC)
  stat_list <- c("AIC", "BIC", "eBIC", "mAIC", "mBIC")
  names(f_list) <- stat_list
  cv_list <- c("outLL")

  stat_sel <- stat_list[na.omit(match(crit_list, stat_list))]
  f_stat <- f_list[stat_sel]
  cv_sel <- cv_list[na.omit(match(crit_list, cv_list))]

  result_stat <- c()
  if(length(stat_sel)){
    for(i in 1:length(stat_sel)){
      result_stat[i] <- f_stat[[i]](NEGMoE_obj)
    }
  }
  names(result_stat) <- stat_sel

  result_cv_m <- c()
  if(length(cv_sel)){
    result_cv <- .NEGMoE_cv(NEGMoE_obj, crit_list = cv_sel, ...)
    result_cv_m <- mean(result_cv)
  }

  result <- matrix(c(result_stat, result_cv_m), nrow = 1)
  colnames(result) <- c(stat_sel, cv_sel)

  return(result)
}

.cvNEGMoE <- function(NEGMoE_obj, g1 = 10, itmax_lambda = 3,
                     itmax_fit = 50, itmax_cv = 20,
                     crit_list = c("all"),
                     lambda1_M = NULL, ...){


  p <- ncol(NEGMoE_obj@X)
  q <- ncol(NEGMoE_obj@Z)
  n <- nrow(NEGMoE_obj@Z)
  K <- NEGMoE_obj@K

  lambda2 = NEGMoE_obj@params$lambda2

  if(is.null(lambda1_M)){
    lambda1_M <- .lambda1Generator(NEGMoE_obj, lambda2, itmax_lambda)
  }

  lambda1_m <- log(p)/(K*n)

  gamma_init <- NEGMoE_obj@NEGMoE_output$gamma
  result = list()

  if("all" %in% crit_list){
    crit_list <- c("AIC", "BIC", "eBIC", "mAIC", "mBIC", "outLL")
  }


  param_temp <- matrix(0, nrow = g1, ncol = 2)
  result_temp <- matrix(0, nrow = g1, ncol = length(crit_list))
  C <- (log(lambda1_M) - log(lambda1_m)) / (g1 - 1)
  lambda1_seq <- lambda1_m * exp((seq(1,g1) - 1)*C)
  NEGMoE_temp <- NEGMoE_obj
  NEGMoE_temp@X <- NEGMoE_temp@X
  NEGMoE_temp@params$verbose = FALSE
  NEGMoE_temp@params$itmax = itmax_fit

  for(j in 1:g1){
    message("Evaluate on lambda1 = ", lambda1_seq[j],"...")
    NEGMoE_temp@params$lambda1 = lambda1_seq[j]
    NEGMoE_temp <- fitNEGMoE(NEGMoE_temp,
                            beta_init = NULL,
                            gamma_init = gamma_init)
    param1 = paste(NEGMoE_temp@params$lambda1, collapse = ",")
    param2 = paste(NEGMoE_temp@params$lambda2, collapse = ",")
    param_temp[j,] <- c(param1, param2)
    res_temp <- calcCriterion(NEGMoE_temp, crit_list = crit_list,
                              itmax = itmax_cv)
    result_temp[j,] <- res_temp
  }
  param_df <- as.data.frame(param_temp)
  colnames(param_df) <- c("lambda1", "lambda2")
  result_df <- as.data.frame(result_temp)
  colnames(result_df) <- colnames(res_temp)
  result <- cbind(param_df, result_df)

  return(result)
}

.chooselambda1 <- function(result, crit_sel = "outLL"){

  lambda1_sel <- c()
  idx_sel <- which.max(result[,crit_sel])
  lambda1_sel <- as.numeric(result[idx_sel,1])

  return(lambda1_sel)
}

###########################################################################
# function to Cross validate NEGMoE with different tunning parameters lambda
##########################################################################
#' Parameters tunning in NEGMoE
#' @description This function calculate evaluation of NEGMoE object.
#' @param NEGMoE_obj a NEGMoE of object without cross validated parameters.
#' @param verbose  A logical input indicating whether the intermediate
#' steps will be printed.
#' @param ... Other parameters that can be passed to cvNEGMoE.
#' @return A NEGMoE object with different lambdas.
#' @export
#' @seealso \code{\link{createCVList}}

cvNEGMoE <- function(NEGMoE_obj, verbose = T, ...){

  lambda2_seq <- NEGMoE_obj@cvParams$lambda2
  g1 <- NEGMoE_obj@cvParams$g1
  itmax_cv <- NEGMoE_obj@cvParams$itmax_cv
  itmax_fit <- NEGMoE_obj@cvParams$itmax_fit
  itmax_lambda <- NEGMoE_obj@cvParams$itmax_lambda
  crit_eval <-  NEGMoE_obj@cvParams$crit_eval
  crit_sel <- NEGMoE_obj@cvParams$crit_sel
  track <- NEGMoE_obj@cvParams$track
  gamma_init <- NEGMoE_obj@NEGMoE_output$gamma

  if("all" %in% crit_eval){
    crit_eval <- c("AIC", "BIC", "eBIC", "mAIC", "mBIC", "outLL")
  }

  result_lambda1 <- list()
  lambda1_sel <- c()
  result_lambda2 <- matrix(0, nrow = length(lambda2_seq),
                           ncol = length(crit_eval))
  param_temp <- matrix(0, nrow = length(lambda2_seq), ncol = 2)

  for(i in 1:length(lambda2_seq)){
    message("Evaluate on lambda2 = ", lambda2_seq[i], ".....")
    NEGMoE_temp <- NEGMoE_obj
    NEGMoE_temp@params$lambda2 = lambda2_seq[i]
    result_lambda1[[i]] <- .cvNEGMoE(NEGMoE_temp, g1, itmax_lambda,
                                    itmax_fit, itmax_cv,
                                    crit_eval)

    lambda1_sel[i] <- .chooselambda1(result_lambda1[[i]], crit_sel = crit_sel)

    NEGMoE_temp@params$lambda1 = lambda1_sel[i]
    NEGMoE_temp@params$verbose = FALSE
    NEGMoE_temp@params$itmax = itmax_fit

    NEGMoE_temp <- fitNEGMoE(NEGMoE_temp,
                            beta_init = NULL,
                            gamma_init = gamma_init)

    result_temp <- calcCriterion(NEGMoE_obj = NEGMoE_temp, crit_list = crit_eval,
                                 itmax = itmax_cv, ...)

    result_lambda2[i,] <- result_temp

    param1 = paste(NEGMoE_temp@params$lambda1, collapse = ",")
    param2 = paste(NEGMoE_temp@params$lambda2, collapse = ",")
    param_temp[i,] <- c(param1, param2)
  }
  param_df <- as.data.frame(param_temp)
  colnames(param_df) <- c("lambda1", "lambda2")
  result_df <- as.data.frame(result_lambda2)
  colnames(result_df) <- colnames(result_temp)
  result_lambda2 <- cbind(param_df, result_df)

  idx_sel <- which.max(result_lambda2[,crit_sel])
  lambda1_choose <- lambda1_sel[[idx_sel]]
  lambda2_choose <- lambda2_seq[idx_sel]

  if(track){

    cv_result <- list(result_lambda1 = result_lambda1,
                      result_lambda2 = result_lambda2,
                      lambda1_choose = lambda1_choose,
                      lambda2_choose = lambda2_choose)
  }else{
    cv_result <- list(result_lambda2 = result_lambda2,
                      lambda1_choose = lambda1_choose,
                      lambda2_choose = lambda2_choose)
  }
  NEGMoE_obj@cvResult <- cv_result

  return(NEGMoE_obj)
}
