#' A wrapper function for Elastic net for each variable (pseudo-likelihood) 
#' @param X data matrix for calculate pseudo-likelihood using elastic net (n x p).
#' @param lambda parameter lambda for elastic net.
#' @param alpha parameter lambda for elastic net(if alpha = 1 for lasso).
#' @param R weight matrix for each observations (n x 1).
#' @param beta_init initial value for beta.
#' @param itmax maximum number of iterations.
#' @param eps A numeric threshold for coordinate descent.
#' @param eps_var A numeric threshold to exclude variables with variance smaller than the threshold
#' @param intercept A logic variable to indicate whether include intercept terms.
#' @param version If version = "glmnet" use glmnet package for fitting or use
#' self-written coordinate descent.
#' @return a matrix of regression coefficients.
#' @export
Elnet_ps <- function(X, lambda, alpha = 1, R = NULL, beta_init = NULL,
                     itmax = 1e3, eps = 1e-3, eps_var = 1e-5,
                     intercept = TRUE, version = "glmnet"){
  
  
  n <- nrow(X)
  p <- ncol(X)
  
  B <- matrix(0, nrow = (p + 1), ncol = p)
  
  if(p == 2){

    fit1 <- lm(X[,1]~X[,2])
    fit2 <- lm(X[,2]~X[,1])
    
    B[-2,1] <- as.matrix(coef(fit1))
    B[-3,2] <- as.matrix(coef(fit2))
    
  }else{
    
    lambda <- ext_param(lambda, p)
    
    for(i in 1:p){
      
      X_temp <- X[,-i]
      Y_temp <- X[,i]
      
      if(version == "glmnet"){
        B[-(i + 1), i] <- .Elnet_glmnet(X = X_temp, y = Y_temp, lambda = lambda[i], alpha = alpha, R = R, beta_init = beta_init,
                               itmax = itmax, eps = eps, eps_var = eps_var, intercept = intercept)
      }else{
        B[-(i + 1), i] <- .Elnet(X = X_temp, y = Y_temp, lambda = lambda[i], alpha = alpha, R = R, beta_init = beta_init,
                               itmax = itmax, eps = eps, eps_var = eps_var, intercept = intercept)
      }
  
      
      
    }
  }
  
  return(B)
  
}