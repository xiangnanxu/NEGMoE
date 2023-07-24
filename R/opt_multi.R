
#' Inner loop (Elastic net) for sparse multinomial regression with proximal Newtown updating.

.sMulti_single <- function(X, y, beta0, residual_k, R, backtracking = T,
                          lambda, alpha, itmax = 1e3, beta_max = 10, eps = 1e-5,
                          version = "my", intercept = TRUE){

  n <- nrow(X)
  p <- ncol(X)

  beta0 <- as.matrix(beta0)
  beta <- beta0

  X_1 <- cbind(rep(1,nrow(X)), X)

  it = 0

  repeat{

    pi_tid <- 1 / (1 + residual_k * exp(- X_1 %*% beta))
    pi_tid[is.nan(pi_tid)] <- 0
    pi_tid[pi_tid < 1e-8] <- 0

    w_tid <- pi_tid * (1 - pi_tid)

    w_tid[w_tid < 1e-8] <- 1e-8

    weight <- as.vector(w_tid * R)

    z <- X_1 %*% beta + (y - pi_tid) / w_tid

    if(version == "glmnet"){

      if(var(z) < 1e-8){

        beta <- matrix(0, nrow = p + 1, ncol = 1)

        beta[1] <- mean(z)
      }else{

        if(p == 1){
          beta <- .Elnet(X, z, lambda, alpha, beta_init = beta0, R = weight, intercept = intercept, itmax = itmax)
        }else{

          beta <- try(.Elnet_glmnet(X, z, lambda, alpha, R = weight, intercept = intercept))

          if("try-error" %in% class(beta)){
            beta <- .Elnet(X, z, lambda, alpha, beta_init = beta0, R = weight, intercept = intercept, itmax = itmax)
          }

        }

        beta <- clipping(beta, beta_max, -beta_max)
      }

    }else{
      beta <- .Elnet(X, z, lambda, alpha, beta_init = beta0, R = weight, intercept = intercept, itmax = itmax)
    }

    if(backtracking){
      obj_logistic0 <- purrr::partial(.obj_logistic, R = R, X = X, y = y,
                                      alpha = alpha, lambda = lambda,
                                      intercept = T, multiply = TRUE)


      beta <- clipping(beta, beta_max, -beta_max)

      beta <- .backtracking_f(obj_logistic0, x_old = beta0, x_new = beta, mm = "min")
    }

    beta <- hard_thresh(beta)

    cond1 <- (sum(abs(beta - beta0)) < eps)
    cond2 <- (it > itmax)

    if(cond1|cond2){
      break()
    }

    it <- it + 1
    beta0 <- beta
  }

  return(beta0)

}

#' Inner loop (Coordinate descent) for sparse multinomial regression with proximal Newtown updating.

.sMulti_step <- function(X, y, beta0, lambda, R = NULL, backtracking = T,
                        alpha, itmax = 1e3, beta_max = 10, version = version){

  n <- nrow(X)
  p <- ncol(X)

  K <- ncol(y)

  if(!length(R)){
    R = rep(1, n)
  }else{
    R = as.vector(R)
  }

  beta <- beta0

  X_1 <- cbind(rep(1,nrow(X)), X)

  for(k in 2:K){

    beta_k <- beta[,-k]
    residual_k <- rowSums(exp(X_1 %*% beta_k))

    beta[,k] <- .sMulti_single(X, y = y[,k], beta0 = beta0[,k], residual_k, R, lambda, alpha,
                              itmax = itmax, beta_max = beta_max,
                              version = version, backtracking = backtracking)

  }

  beta <- clipping(beta, beta_max, -beta_max)

  return(beta)
}

# Outer loop for sparse multinomial regression

.sMulti1 <- function(X, y ,lambda, alpha = 1, R = NULL, beta_init = NULL,
                   itmax1 = 1e3, itmax2 = 10, eps = 1e-8, backtracking = T,
                   beta_max = 10, intercept = TRUE, version = "my"){

  y <- as.matrix(y)

  n <- nrow(X)
  p <- ncol(X)

  K <- ncol(y)

  lambda <- lambda

  if(K == 1){
    K = 2
    y = cbind(1 - y, y)
  }

  if(!length(R)){
    R = rep(1, n)
  }else{
    R = as.vector(R)
  }

  if(!length(beta_init)){

    beta <- matrix(0, nrow = (p + 1), ncol = K)
    for(i in 2:K){

      temp <- sum(y[,i] == 1)

      if(temp < 7){
        beta[1,i] = -beta_max
      }else if(temp > (n - 7)){
        beta[1,i] = beta_max
      }else{
        temp <- temp/n

        beta[1,i] = log(temp/(1 - temp))
      }
    }

  }else{
    if(K == 2){
      beta <- matrix(0, nrow = (p + 1), ncol = K)
      if(ncol(beta_init) == 1){
        beta[,2] <- beta_init
      }
      else{
        beta_init <- V_trans(beta_init)
        beta <- beta_init[,1:2]
      }
    }else{
      beta_init <- V_trans(beta_init)
      beta <- beta_init
    }
  }

  X_1 <- cbind(rep(1,nrow(X)), X)

  beta_old <- beta

  it = 0

  obj_smulti <- c()

  repeat{

    beta <- .sMulti_step(X, y, beta, R = R, lambda = lambda, alpha = alpha,
                        itmax = itmax2, beta_max = beta_max,
                        version = version, backtracking = backtracking)

    cond1 <- sum(abs(beta - beta_old)) < eps

    cond2 <- (it > itmax1)
    if(cond1|cond2){

      break()

    }

    it = it + 1
    beta_old <- beta

  }


  return(beta)

}

#' Wrapper function for sparse Multinomial regression
#' @param X description
#' @param y description
#' @param R description
#' @param lambda description
#' @param alpha description
#' @param beta_max description
#' @param intercept description
#' @param var_eps description
#' @return description

sMulti_glmnet <- function(X, y, R = NULL,lambda, alpha = 1, beta_max = 10,
                          intercept = TRUE, var_eps = 1e-5){

  X <- as.matrix(X)

  n <- nrow(X)
  p <- ncol(X)

  y <- as.matrix(y)

  y[is.infinite(y)] <- 0
  y[is.na(y)] <- 0

  K <- ncol(y)

  if(!length(R)){
    R = rep(1, n)
  }else{
    R = as.vector(R)
  }

  lambda = lambda

  beta <- matrix(0, nrow = (p + 1), ncol = K)

  lambda <- lambda

  if(p == 1){

    beta <- sMulti(X = X, y = y, lambda = lambda, alpha = alpha,  R = R, beta_max = beta_max, intercept = T, version = "my")

  }else{

    if(K == 1){

      if(!(class_check(y))){
        temp <- mean(y)
        temp_beta0 <- log(temp/(1 - temp))
        temp_beta0 <- clipping(temp_beta0, beta_max, -beta_max)
        beta[1,1] <- temp_beta0
        beta[2:(p + 1),1] <- 0

      }else{
        beta <- suppressWarnings(as.matrix(coef(glmnet::glmnet(X, y, family = "binomial", weights = R,
                                                               lambda = lambda * n / sum(R),
                                                               alpha = alpha, standardize = F))))
        beta <- clipping(beta, beta_max, -beta_max)

      }
    }else{

      if(!(class_check(y))){

        beta <- .sMulti(X = X, y = y, lambda = lambda, alpha = alpha,  R = R, beta_max = beta_max, intercept = T, version = "my")

      }else{

        glmtemp <- try(glmnet::glmnet(X, y, family = "multinomial", lambda = lambda * n / sum(R),
                                      weights = R, alpha = alpha, standardize = F))

        if("try-error" %in% class(glmtemp)){

          beta <- .sMulti(X = X, y = y, lambda = lambda, alpha = alpha,  R = R, beta_max = beta_max,
                         intercept = T, version = "my")

        }else{

          beta_temp <- suppressWarnings(coef(glmtemp))

          for(i in 1:K){
            beta[,i] <- as.matrix(beta_temp[[i]])
          }

          beta <- V_trans(beta,beta_max = beta_max)
        }



      }
    }
  }
  beta <- as.matrix(beta)
  return(beta)
}

.sMulti <- function (X, y, lambda, alpha = 1, R = NULL, beta_init = NULL, 
                    itmax1 = 1000, itmax2 = 10, eps = 1e-08, backtracking = T, 
                    beta_max = 10, intercept = TRUE, version = "my") 
{
  y <- as.matrix(y)
  n <- nrow(X)
  p <- ncol(X)
  K <- ncol(y)
  lambda <- lambda
  if (K == 1) {
    K = 2
    y = cbind(1 - y, y)
  }
  if (!length(R)) {
    R = rep(1, n)
  }
  else {
    R = as.vector(R)
  }
  if (!length(beta_init)) {
    beta <- matrix(0, nrow = (p + 1), ncol = K)
    for (i in 2:K) {
      temp <- sum(y[, i] == 1)
      if (temp < 7) {
        beta[1, i] = -beta_max
      }
      else if (temp > (n - 7)) {
        beta[1, i] = beta_max
      }
      else {
        temp <- temp/n
        beta[1, i] = log(temp/(1 - temp))
      }
    }
  }
  else {
    if (K == 2) {
      beta <- matrix(0, nrow = (p + 1), ncol = K)
      if (ncol(beta_init) == 1) {
        beta[, 2] <- beta_init
      }
      else {
        beta_init <- V_trans(beta_init)
        beta <- beta_init[, 1:2]
      }
    }
    else {
      beta_init <- V_trans(beta_init)
      beta <- beta_init
    }
  }
  X_1 <- cbind(rep(1, nrow(X)), X)
  beta_old <- beta
  it = 0
  obj_smulti <- c()
  repeat {
    beta <- .sMulti_step(X, y, beta, R = R, lambda = lambda, 
                        alpha = alpha, itmax = itmax2, beta_max = beta_max, 
                        version = version, backtracking = backtracking)
    cond1 <- sum(abs(beta - beta_old)) < eps
    cond2 <- (it > itmax1)
    if (cond1 | cond2) {
      (break)()
    }
    it = it + 1
    beta_old <- beta
  }
  return(beta)
}

