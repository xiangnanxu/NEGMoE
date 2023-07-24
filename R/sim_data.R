#' Generate simulated data
#' @param n Number of samples
#' @param p Number of variables of data matrix X (input of experts network)
#' @param q Number of variables of data matrix Z (input of gating network)
#' @param s Sparsity parameters, control non-zero coefficients in gating network coefficients.
#' @param K Number of mixture components
#' @param niu_beta Strength of beta for gating network
#' @param r Separation parameters for latent class. Larger r encourage more separation in latent classes.
#' @param inv_conn Number of non-zero components in inverse covariance matrix.
#' @param version link function, could be "logit" or "probit"
#' @return A list contain simulated data, including data matrix, X, input of gating network Z, coefficients of 
#' gating network beta_V, real mean of X mu, and real inverse covariance matrix of X Theta.
dat_gen <- function(n, p = 20, q = 20, s = 5, K = 2, niu_beta = 2, r= 0.1,
                    inv_signal = 0.8, inv_conn = 10, version = "probit"){
  
  Z <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = diag(rep(1, q)))
  
  q_in1 <- integer(0)
  q_out1 <- 1:q
  
  beta <- matrix(0, nrow = q, ncol = K)
  
  for(i in 1:K){
    idx_gen = .data_gen_idx(q_in1, q_out1, s)
    
    q_in1 <- idx_gen$p_in
    q_out1 <- idx_gen$p_out
    
    beta[,i] = .data_gen_beta(q, idx_gen$idx, niu = niu_beta)
  }
  
  
  pi_temp <- Z %*% beta
  
  pi <- .data_gen_latent(Z, K, beta, r)
  
  if(version == "logit"){
    latent <- sample_prob(pi)
  }else{
    latent <- hard_prob(pi)
  }
  
  latent <- hard_prob(pi)
  
  Theta <- array(0, dim = c(p, p, K))
  S <- array(0, dim = c(p, p, K))
  
  sel_list <- combn(p, 2)
  mu <- matrix(0, nrow = K, ncol = p)
  for(i in 1:K){
    
    comb_sel <- sample(1:ncol(sel_list), inv_conn, replace = F)
    
    Theta_temp <- diag(rep(1,p))
    
    for(j in 1:inv_conn){
      
      val_temp <- rnorm(1, mean = 0.6, sd = 0.1) * sign(runif(1, min = -1, max = 1))
      
      Theta_temp[sel_list[1,comb_sel[j]], sel_list[2,comb_sel[j]]] <- val_temp
      Theta_temp[sel_list[2,comb_sel[j]], sel_list[1,comb_sel[j]]] <- val_temp
      
    }
    
    Theta_temp <- proj_sd(Theta_temp, thresh = 0.2)
    
    
    Theta[,,i] <- Theta_temp
    S[,,i] <- solve(Theta_temp)
    
    mu_temp <- rnorm(p)
    mu[i,] <- mu_temp/sqrt(sum(mu_temp^2))
    
  }
  
  
  X <- matrix(0, nrow = n, ncol = p)
  
  for(i in 1:K){
    
    idx <- (latent[,i] == 1)
    sub_n <- sum(idx)
    
    X[latent[,i] == 1,] = MASS::mvrnorm(n = sub_n, mu = mu[i,], Sigma = S[,,i])
  }
  
  
  return(list(X= X, Z= Z, latent = latent, pi = pi, S = S, Theta = Theta, mu = mu, beta_V = beta))
}

.data_gen_beta <- function (p, idx, niu = 0.3) 
{
  beta <- rep(0, p)
  beta_temp <- rnorm(idx, sd = 1)
  beta[idx] <- beta_temp + sign(beta_temp) * niu
  return(beta)
}

.data_gen_latent <- function (Z, K, beta_V, r = 0.2) 
{
  n <- nrow(Z)
  p1 <- ncol(Z)
  latent_odd <- Z %*% beta_V + r * MASS::mvrnorm(n, mu = rep(0, 
                                                             K), 
                                                 Sigma = diag(rep(1, K))) * sqrt(sum(abs(beta_V)^2))
  latent_prob <- t(apply(latent_odd, MARGIN = 1, softmax))
  return(latent_prob)
}

.data_gen_idx <- function (p_in, p_out, s) 
{
  p <- length(p_in) + length(p_out)
  if (s > 1) {
    s_real <- s
  }
  else {
    s_real <- ceiling(p * s)
  }
  if (!length(p_out)) {
    idx <- sample(p_in, s_real)
  }
  else if (length(p_out) < s_real) {
    idx <- sample(p_in, s_real - length(p_out))
    idx <- union(idx, p_out)
  }
  else {
    idx <- sample(p_out, s_real)
  }
  p_in <- union(p_in, idx)
  p_out <- setdiff(p_out, idx)
  return(list(p_in = p_in, p_out = p_out, idx = idx))
}