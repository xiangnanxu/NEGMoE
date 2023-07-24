#' Graphical lasso algorithm
#' @param X A data matrix of n observations and p variables.
#' @param mu A vector of means for each variable(p x 1).
#' @param lambda penalty parameters for Glasso.
#' @param weight weighted vector for each observations (n x 1). If is NULL, will use all 1 vectors.
#' @param itmax Maximum number of iterations. By default is 1e2.
#' @param precise A logic variable indicate whether to use adaptive or ordinary graphical lasso.
#' @param Theta_init Initial value of Theta. If is NULL, will use a diagonal matrix.
#' @return An estimation of inverse covariance matrix (Theta).
GLasso <- function(X, mu = NULL, lambda, weight = NULL,
                   itmax = 1e2, precise = T, Theta_init = NULL){
  
  n <- nrow(X)
  p <- ncol(X)
  
  if(is.null(mu)){
    mu = colMeans(X)
  }
  
  
  if(is.null(weight)){
    weight <- rep(1, n)
  }
  
  if(sum(weight) < 1e-8){
    return(diag(rep(1, p)))
  }
  
  S <-  (t(X - mu) %*% diag(weight) %*% (X - mu))/sum(weight)
  
  W = S
  beta = -lambda * diag(rep(1,p))
  if(is.null(Theta_init)){
    Theta = diag(rep(1,p))
  }else{
    Theta = Theta_init
  }

  it = 1
  
  LL = c(.LLglasso(X, mu, Theta, weight = weight))
  PLL = c(.PLLglasso(X, mu, Theta, weight = weight, lambda = lambda))
  cond_it = c(F)
  cond_PLL = c(F)
  cond_dev = c(F)
  
  repeat{
    
    for (i in 1:p){
      
      idx <- 1:p
      
      idx <- idx[-i]
      
      beta[idx, i] <- .myglasso(W[idx,idx], S[idx,i], lambda, beta[idx,i])
      
      W[idx,i] <- W[idx, idx] %*% beta[idx,i]
      
    }

    
    it = it + 1
    
    for(j in 1:p){
      
      idx <- 1:p
      
      idx <- idx[-j]
      
      Theta[j,j] <- 1/(W[j,j] - sum(W[idx,j] * beta[idx,j]))
      Theta[idx, j] <- (-Theta[j,j] * beta[idx,j]) 
      
    }
    
    LL[it] <- .LLglasso(X, mu, Theta, weight = weight)
    PLL[it] <- .PLLglasso(X, mu, Theta, weight = weight, lambda = lambda)
    
    cond_it[it] <- .stop_it(it, itmax = itmax)
    
    cond_PLL[it] <- .stop_LL(PLL, eps = 1e-2)
    
    cond_dev[it] <- .stop_dev(W, S, lambda, beta, eps = 1e-5)
    
    if(precise){
      if((cond_it[it])|(cond_PLL[it])){
        break()
      }
    }else{
      if((cond_it[it])|(cond_PLL[it])|(cond_dev[it])){
        break()
      }
    }
    
    
  }
  
  cond <- data.frame(it = cond_it, LL = cond_PLL, KKT = cond_dev)
  return(list(Theta= Theta, LL = LL, PLL= PLL, stop_cond = cond))
  
}

#' Calculate Likelihood function for gglasso
.LLglasso <- function(X, mu = NULL, Theta, weight = NULL, reduce = T, type = "icov"){
  
  if(type == "cov"){
    Theta = chol_solve(Theta)
  }
  
  if(is.null(mu)){
    mu = colMeans(X)
  }
  
  n <- nrow(X)
  p <- ncol(X)
  
  if(is.null(weight)){
    weight <- rep(1, n)
  }
  
  if(sum(weight) < 1e-8){
    return(diag(rep(1, p)))
  }
  
  det_term <- log(det(Theta)) * weight
  
  X_temp <- (X - mu) * sqrt(weight)
  
  trace_term <- apply(X_temp, 1, function(vec){
    A = outer(vec, vec, "*")
    B = A * Theta
    return(sum(B))
  })
  
  if(reduce){
    return(sum(det_term - trace_term))
  }else{
    return(det_term - trace_term)
  }
  
}

#' Calculate penalized log-likelihood for graphical lasso.
.PLLglasso <- function(X, mu = NULL, Theta, lambda, weight = NULL, type = "icov"){
  
  n <- nrow(X)
  p <- ncol(X)
  
  if(type == "cov"){
    Theta = chol_solve(Theta)
  }
  
  if(is.null(mu)){
    mu = colMeans(X)
  }
  
  if(is.null(weight)){
    weight <- rep(1, n)
  }
  
  if(sum(weight) < 1e-8){
    return(diag(rep(1, p)))
  }
  
  n <- nrow(X)
  p <- ncol(X)
  

  LL <- .LLglasso(X = X, mu = mu, Theta = Theta, weight = weight, reduce = T)
  
  PLL <- LL - lambda * sum(abs(Theta[upper.tri(Theta)])) * sum(weight)

  return(PLL)
  
}

# Stopping criteria based on log-likelihood
.stop_LL <- function(LL, eps = 1e-9){
  l <- length(LL)
  if(l > 1){
    
    return(abs(LL[l - 1] - LL[l]) < eps)
    
  }else{
    
    return(FALSE)
    
  }
}

# Stopping criteria based on each variables.

.stop_dev <- function(W, S, lambda, beta, eps = 1e-9){
  
  p <- nrow(W)
  
  cond <- matrix(F, nrow = p, ncol = p)
  
  idx <- (abs(beta) < 1e-9)
  
  cond[idx] <- abs(W[idx] - S[idx]) <= lambda
  cond[!idx] <- (abs(W[!idx] - S[!idx] + lambda * sign(beta[!idx])) < eps)
  diag(cond) = T
  
  
  return(all(cond))
}

# Self-written gglaso algorithm for gglasso
.myglasso <- function(W, S, lambda, beta, eps = 1e-9, itmax = 1e5){
  
  p <- ncol(W)
  
  df <- W %*% beta - S
  
  cond <- rep(F, p)
  
  it = 0
  
  repeat{
    
    for(j in 1:p){
      
      b <- beta[j]
      
      beta[j] <- soft_thresh(b - df[j] / W[j,j], lambda / W[j,j])
      
      if(beta[j] != b){
        df <- df + (beta[j] - b) * W[,j]
      }
      
    }
    
    it = it + 1
    
    if(it > itmax){
      break()
    }else{
      
      idx <- (beta == 0)
      cond[idx] <- (abs(df[idx]) <= lambda)
      cond[!idx] <- (abs(df[!idx] + lambda * sign(beta[!idx])) < eps)
      if(all(cond)){
        break()
      }
    }
    
  }
  return(beta)
  
}

#' soft threshold function
#' @param x input of vector for soft-threshold function
#' @param lambda penalty parameters for soft-threshold function.
#' @return output is sign(x)*max(x-lambda,0)
soft_thresh <- function(x, lambda){
  
  soft_x <- sign(x) * max(0, abs(x) - lambda)
  
  return(soft_x)
}

#' All kinds of likelihood function for graphical lasso.
#' @param X Input matrix of X for experts network.
#' @param Z Input matrix of Z for gating network.
#' @param beta_V Coefficients matrix for gating network.
#' @param Theta An array for estimated inverse covariance matrices.
#' @param mu mean vector for each experts.
#' @param r_i likelihood ratio for each components.
#' @param lambda1 Penalties for graphical lasso.
#' @param lambda2 Penalties for gating networks.
#' @param ratio_max clipping for the maximium of ratio to be in one cluster.
#' @return A vector contain Log-likelihood, Penalized Log-likelihood, Complete log-likelihood and Penalized complete log-likelihood.
LL_Glasso <- function(X, Z, beta_V, Theta, mu, r_i, lambda1, lambda2, ratio_max = 20){
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  
  K <- ncol(beta_V)
  X_1 <- cbind(rep(1,nrow(X)), X)
  Z_1 <- cbind(rep(1,nrow(Z)), Z)
  
  LL <- c()
  PLL <- c()
  
  LL_comp <- c()
  PLL_comp <- c()
  
  pi <- prob_multi(Z_1, beta_V)
  
  Likeli_G <- matrix(0, nrow = n, ncol = K)
  
  penalty1 <- c()
  
  for(i in 1:K){
    
    Theta_temp <- Theta[,,i]
    Likeli_G[,i] <- .LLglasso(X = X, mu = mu[i,], Theta = Theta_temp, reduce = F)
    penalty1[i] <- lambda1[i] * sum(abs(Theta_temp[upper.tri(Theta_temp)])) * sum(r_i[,i])
    
  }
  
  penalty2 <- lambda2 * sum(abs(beta_V[2:nrow(beta_V),])) * n
  
  pi_ratio <- log(pi) - log(pi[,1])
  Likeli_ratio <- Likeli_G - Likeli_G[,1]
  
  total_ratio <- pi_ratio + Likeli_ratio
  
  total_ratio[total_ratio > ratio_max] = ratio_max
  total_ratio[total_ratio < -ratio_max] = -ratio_max
  
  LL <- Likeli_G[,1] + log(pi[,1]) + log(rowSums(exp(total_ratio)))
  LL0 <- sum(LL)
  PLL0 <- LL0 - sum(penalty1) - penalty2
  
  LL_comp <- sum(r_i * Likeli_G)
  PLL_comp <- LL_comp - sum(penalty1) - penalty2
  
  return(c(LL0, PLL0, LL_comp, PLL_comp))
}

#' Project a matrix to semi-definite space.
#' @param X Input of the matrix
#' @param thresh threshold for eigenvalue of output semidefinite matrix.
#' @return A semi-definite matrix by iterative projection.
#' @export
proj_sd <- function(X ,thresh = 0.1){
  
  X_0 <- (X + t(X))/2
  it = 0
  
  repeat{
    
    eig_X <- eigen(X_0)
    
    U <- eig_X$vectors
    lambda <- eig_X$values
    
    pmax <- purrr::partial(max, thresh)
    lambda <- sapply(lambda, pmax)
    
    X_1 <- U %*% diag(lambda) %*% t(U)
    
    X_1[X == 0] = 0
    
    
    cond_conv <- ((sum(abs(X_1 - X_0))/sum(abs(X_0))) < 1e-8)
    cond_it <- .stop_it(it, itmax = 1e5)
    
    if((cond_conv)|(cond_it)){
      break()
    }
    
    it = it + 1
    
    X_0 <- X_1
  }
  return(X_1)
}

#' Project a matrix to correlation space with non-zero correlation larger than thresh_cor.
#' @param X Input of the matrix
#' @param thresh_cor threshold for minimal non-zero correlation of output semi-definite matrix.
#' @param thresh_eig threshold for eigenvalue of output semi-definite matrix.
#' @return A semi-definite matrix by iterative projection.
proj_sd_thresh <- function(X, thresh_cor = 0.3, thresh_eig = 0.1){
  
  X_0 <- (X + t(X))/2
  it = 0
  
  thresh_inv_eig <- 1/thresh_eig
  
  repeat{
    
    R_0 <- diag(1/sqrt(diag(X_0))) %*% X_0 %*% diag(1/sqrt(diag(X_0)))
    
    X_0[abs(R_0) < thresh_cor] = 0
    
    eig_X <- eigen(X_0)
    U <- eig_X$vectors
    lambda <- eig_X$values
    
    lambda[lambda < thresh_eig] <- thresh_eig
    lambda[lambda > thresh_inv_eig] <- thresh_inv_eig
    
    X_1 <- U %*% diag(lambda) %*% t(U)

    cond_conv <- ((sum(abs(X_1 - X_0))/sum(abs(X_0))) < 1e-8)
    cond_it <- .stop_it(it, itmax = 1e5)
    
    if((cond_conv)|(cond_it)){
      break()
    }
    
    it = it + 1
    
    X_0 <- X_1
    
  }
  
  return(X_0)  
}

#' Inverse matrix solver use Cholesky decomposition
#' @param S A semi-definite matrix.
#' @return S^{-1}
chol_solve <- function(S){
  
  L <- chol(S)
  
  L_inv <- solve(L)
  
  return(L_inv %*% t(L_inv))
  
}
