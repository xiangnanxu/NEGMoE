#' A wrapper function for Elastic net for multivariate response (Y multivariate)
#' @param X covariate for elastic net regression (n x p).
#' @param Y response variable for elastic net regression (n x q).
#' @param lambda parameter lambda for elastic net.
#' @param alpha parameter lambda for elastic net(if alpha = 1 for lasso).
#' @param R weight matrix for each observations (n x 1).
#' @param beta_init initial value for beta.
#' @param itmax maximum number of iterations. By default is 1e3.
#' @param eps A numeric threshold for coordinate descent. By default is 1e-3.
#' @param eps_var A numeric threshold to exclude variables with variance smaller than the threshold.
#' By default is 1e-5.
#' @param intercept A logic variable to indicate whether include intercept terms.
#' @param version If version = "glmnet" use glmnet package for fitting or use
#' self-written coordinate descent.
#' @return a matrix of regression coefficients.          

Elnet_multi <- function(X, Y, lambda, alpha = 1, R = NULL, beta_init = NULL,
                       itmax = 1e3, eps = 1e-3, eps_var = 1e-5, intercept = TRUE,
                       version = "glmnet"){

  X <- as.matrix(X)
  Y <- as.matrix(Y)

  n <- nrow(X)
  p <- ncol(X)

  nr <- ncol(Y)

  X_1 <- cbind(rep(1, n), X)

  lambda <- ext_param(lambda, nr)

  if(is.null(R)){
    R <- rep(1, n)
  }

  if(p == 1){

    fit_lm <- lm(Y~X, weights = R)
    B <- as.matrix(coef(fit_lm))

    return(B)
  }else{

    if(version == "glmnet"){

      B <- matrix(0, nrow = (p + 1), ncol = nr)

      for(i in 1:nr){

        B[,i] <- .Elnet_glmnet(X = X, y = Y[,i], lambda = lambda[i], alpha = alpha,
                              R = R, intercept = intercept, beta_init = beta_init,
                              eps = eps, eps_var = eps_var, itmax = itmax)
      }

    }else{

      B <- matrix(0, nrow = (p + 1), ncol = nr)

      for(i in 1:nr){

        B[,i] <- .Elnet(X = X, y = Y[,i], lambda = lambda[i], alpha = alpha,
                              R = R, intercept = intercept, beta_init = beta_init,
                              eps = eps, eps_var = eps_var, itmax = itmax)
      }

    }

  }

  return(B)


}
