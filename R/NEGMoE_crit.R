###################################################################
# function to predict output of new data using fitted NEGMoE object
###################################################################
#' Make predictions from a fitted "NEGMoE" object.
#' @description This function predict new data using.
#' @param NEGMoE_obj A NEMoE_obj object with fitted parameters.
#' @param X_new A matrix of new X.
#' @param Z_new A matrix of new Y.
#' @param Y_new A matrix of new Z.
#' @param full Whether output the result in both gating network
#'  and experts network.
#' @return A number of fitted (pseudo) Likelihood or predicted response.
#' @export
NEGMoE_predict <- function(NEGMoE_obj, X_new, Z_new, Y_new = NULL, out= "likelihood", full = TRUE){
  
  type = NEGMoE_obj@type
  K = NEGMoE_obj@K
  gamma = NEGMoE_obj@NEGMoE_output$gamma
  standardize = NEGMoE_obj@standardize
  
  if(is.vector(X_new)){
    X_new <- matrix(X_new, nrow = 1)
  }
  
  if(is.vector(Z_new)){
    Z_new <- matrix(Z_new, nrow = 1)
  }
  
  if(is.vector(Y_new)){
    Y_new <- matrix(Y_new, nrow = 1)
  }
  
  if(standardize){
    X_new <- (X_new - NEGMoE_obj@.transformation$mu_X) /(NEGMoE_obj@.transformation$sd_X) 
    Z_new <- (Z_new - NEGMoE_obj@.transformation$mu_Z) /(NEGMoE_obj@.transformation$sd_Z) 
  }
  
  n <- nrow(Z_new)
  Z_new1 <- cbind(rep(1,n), Z_new)
  tmp_new <- Z_new1 %*% gamma
  pi_new <- t(apply(tmp_new,1, FUN = softmax))
  
  if(type == "icov"){
  
      experts_Theta = NEGMoE_obj@NEGMoE_output$experts$Theta
      experts_mu = NEGMoE_obj@NEGMoE_output$experts$mu
      
      LL_tmp <- matrix(0, nrow = n, ncol = K)
      
      for(i in 1:K){
        LL_tmp[,i] <- .LLglasso(X_new, mu = experts_mu[i,],
                                Theta = experts_Theta[,,i], reduce = F)
      }
    
  }else if(type == "pseudoLL"){
      experts_W <- NEGMoE_obj@NEGMoE_output$experts$W
      experts_sigma <- NEGMoE_obj@NEGMoE_output$experts$sgm
      LL_tmp <- matrix(0, nrow = n, ncol = K)
      for(i in 1:K){
        LL_tmp[,i] <- .LLikeli_mvr(X_new, experts_W[,,i], Y = X_new,
                                   experts_sigma[i,], intercept = T, reduce = F)
      }
    
  }else{
      experts_W <- NEGMoE_obj@NEGMoE_output$experts$W
      experts_sigma <- NEGMoE_obj@NEGMoE_output$experts$sgm
      LL_tmp <- matrix(0, nrow = n, ncol = K)
      if(out == "likelihood"){
        for(i in 1:K){
          LL_tmp[,i] <- .LLikeli_mvr(X_new, experts_W[,,i], Y = X_new,
                                     experts_sigma[i,], intercept = T, reduce = F)
        }
      }else{
        for(i in 1:K){
          LL_tmp[,i] <- cbind(rep(1,n),X_new) %*% experts_W[,,i]
        }
      }
      
  }
  if(full){
    return(pi_new * LL_tmp)
  }else{
    return(sum(pi_new * LL_tmp))
  }
  
}

###################################################################
# function to calculate AIC of NEGMoE fitting result
###################################################################

.NEGMoE_AIC <- function(NEGMoE_obj, modified = FALSE, LL = "obs"){
  
  K <- NEGMoE_obj@K
  LL_out <- NEGMoE_obj@NEGMoE_output$obj
  
  if(LL == "comp"){
    LogLikeli <- LL_out[nrow(LL_out),3]
  }else{
    LogLikeli <- LL_out[nrow(LL_out),1]
  }
  
  if(modified){
    
    NEGMoE_df = calcdf(NEGMoE_obj, output = "comp")
    df_gamma <- NEGMoE_df$df_gamma
    df_beta <- NEGMoE_df$df_beta
    
    AIC_est <- (-2*LogLikeli) + 2*(df_gamma + sum(df_beta))
    
  }else{
    NEGMoE_df = calcdf(NEGMoE_obj)
    df_gamma <- NEGMoE_df$df_gamma
    df_beta <- NEGMoE_df$df_beta
    
    AIC_est <- (-2*LogLikeli) + 2*(df_gamma + sum(df_beta))
  }
  
  return(-AIC_est)
  
}


.outLL <- function(NEGMoE_obj, X_new, Z_new, Y_new = NULL){
  return(NEGMoE_predict(NEGMoE_obj, X_new, Z_new, Y_new, full = FALSE))
}


###################################################################
# function to calculate BIC of NEGMoE fitting result
###################################################################

.NEGMoE_BIC <- function(NEGMoE_obj, modified = FALSE, type = "BIC", LL = "obs"){
  
  K <- NEGMoE_obj@K
  LL_out <- NEGMoE_obj@NEGMoE_output$obj
  if(LL == "comp"){
    LogLikeli <- LL_out[nrow(LL_out),3]
  }else{
    LogLikeli <- LL_out[nrow(LL_out),1]
  }
  
  n <- nrow(NEGMoE_obj@Z)
  Z1 <- cbind(rep(1,n), NEGMoE_obj@Z)
  tmp <- Z1 %*% NEGMoE_obj@NEGMoE_output$gamma
  pi <- t(apply(tmp, 1, softmax))
  pi1 <- colSums(pi)
  q <- ncol(NEGMoE_obj@Z)
  p <- ncol(NEGMoE_obj@X)
  
  if(modified){
    
    NEGMoE_df = calcdf(NEGMoE_obj, output = "comp")
    df_gamma <- NEGMoE_df$df_gamma
    df_beta <- NEGMoE_df$df_beta
    
    if(type == "eBIC"){
      BIC_est <- (-2 * LogLikeli) + df_gamma*log(n) + df_gamma*log(q)
      for(i in 1:K){
        BIC_est <- BIC_est + df_beta[i]*(log(pi1[i]) + log(p))
      }
      BIC_est <- unname(BIC_est)
    }else{
      df_beta1 <- colSums(df_beta)*log(pi1)
      BIC_est <- (-2 * LogLikeli) + df_gamma*log(n) + sum(df_beta1)
    }
    
  }else{
    
    NEGMoE_df = calcdf(NEGMoE_obj)
    df_gamma <- NEGMoE_df$df_gamma
    df_beta <- NEGMoE_df$df_beta
    
    if(type == "eBIC"){
      BIC_est <- (-2 * LogLikeli) + df_gamma*(log(n) + log(q)) +
        sum(df_beta)*log(n) + sum(df_beta*log(p))
    }else{
      BIC_est <- (-2 * LogLikeli) + (df_gamma + sum(df_beta))*log(n)
    }
  }
  return(-BIC_est)
}

###################################################################
# function to calculate Cross Validation metrics of NEGMoE object
###################################################################

.NEGMoE_cv <- function(NEGMoE_obj, crit_list = "outLL",
                       fold = 5, itmax = 1e2, verbose = F){
  
  p <- ncol(NEGMoE_obj@X)
  q <- ncol(NEGMoE_obj@Z)
  n <- nrow(NEGMoE_obj@Z)
  num_fold <- ceiling(n/fold)
  idx_all <- seq_len(n)
  type = NEGMoE_obj@type
  result <- c()
  
  for(i in 1:fold){
    
    idx_test <- seq((i - 1)*num_fold + 1, min(n, i*num_fold))
    idx_train <- idx_all[-idx_test]
    
    Z_sub = NEGMoE_obj@Z[idx_train,]
    Z_test = NEGMoE_obj@Z[idx_test,]
    X_sub = NEGMoE_obj@X[idx_train,]
    X_test = NEGMoE_obj@X[idx_test,]
    if(!is.null(NEGMoE_obj@Y)){
      Y_sub = NEGMoE_obj@Y[idx_train,]
      Y_test = NEGMoE_obj@Y[idx_test,]
    }else{
      Y_sub = Y_test = NULL
    }
    
    NEGMoE_sub = NEGMoE_obj
    NEGMoE_sub@Z = Z_sub
    NEGMoE_sub@X = X_sub
    NEGMoE_sub@Y = Y_sub
    
    NEGMoE_sub@params$verbose = verbose
    NEGMoE_sub@params$itmax = itmax
    gamma_init = NEGMoE_sub@NEGMoE_output$gamma
    if(type == "icov"){
      beta_init = NEGMoE_sub@NEGMoE_output$experts$Theta
    }else{
      beta_init = NEGMoE_sub@NEGMoE_output$experts$W
    }
    
    NEGMoE_sub = fitNEGMoE(NEGMoE_sub, beta_init = beta_init, gamma_init = gamma_init)
    LL_pred = .outLL(NEGMoE_sub, X_test, Z_test, Y_new = Y_test)
    
    result[i] = LL_pred/length(idx_test)
  }
  return(result)
}
