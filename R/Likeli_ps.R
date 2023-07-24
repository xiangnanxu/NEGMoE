#' Likelihood function for multivariate normal distribution (LL(Y-XB))

.LLikeli_mvr <- function(X, B, Y = NULL, sgmY = NULL, R = NULL, intercept = T, reduce = T){
  
  X <- as.matrix(X)
  
  n <- nrow(X)
  p <- ncol(X)
  B <- as.matrix(B)
  
  if(is.null(Y)){
    Y <- X
  }else{
    Y <- as.matrix(Y)
  }
  
  nr <- ncol(Y)
  
  if(is.null(R)){
    
    R <- rep(1,n)
    
  }
  
  if(is.null(sgmY)){
    
    sgmY <- .est_sgm(X = X, Y = Y, B = B, r_i = R)
    
  }
  
  
  sgm_term <- (- sum(log(sgmY)))
  
  if(intercept){
    
    X <- cbind(rep(1,n), X)

    trace_term <- (-(Y - X %*% B)^2 / (2 * sgmY^2))
    
  }else{
    
    trace_term <- (-(Y - X %*% B)^2 / (2 * sgmY^2))
  }
  
  
  LL <- sgm_term + rowSums(trace_term)
  
  if(reduce){
    
    return(sum( R * LL))  
  }else{
    
    return(R*LL)
  }
}

#' Penalized likelihood function for multivariate normal distribution (PLL(Y-XB))
.PLL_mvr <- function(X, B, Y = NULL, sgmY = NULL, R = NULL, lambda, alpha = 1,
                     intercept = T, multiply = T, div_sgm = T){
  
  X <- as.matrix(X)
  B <- as.matrix(B)
  
  n <- nrow(X)
  p <- ncol(X)
  
  if(is.null(Y)){
    Y <- X
  }else{
    Y <- as.matrix(Y)
  }
  nr <- ncol(Y)
  
  nr <- ncol(Y)
  
  if(is.null(R)){
    
    R <- rep(1,n)
    
  }
  
  if(multiply){
    lambda <- lambda * n
  }
  
  if(length(lambda) == 1){
    lambda <- rep(lambda, nr)
  }
  
  if(is.null(sgmY)){
    
    sgmY <- .est_sgm(X = X, Y = Y, B = B, r_i = R)
    
  }
  
  
  if(div_sgm){
    lambda <- lambda/(sgmY^2)
  }
  
  
  
  LL <- .LLikeli_mvr(X = X, B = B, Y = Y, sgmY = sgmY, R = R, intercept = intercept, reduce = T)
  
  if(intercept){
    
    penalty_temp <- c()
    for(i in 1:nr){
      penalty_temp[i] <-  lambda[i] * (0.5 * (1 - alpha) * sum(B[2:(p + 1), i]^2) + alpha * sum(abs(B[2:(p + 1), i])))
    }
    
    penalty <-  sum(penalty_temp)
    
  }else{
    penalty_temp <- c()
    for(i in 1:nr){
      penalty_temp[i] <- lambda[i] * (0.5 * (1 - alpha) * sum(B[,i]^2) + alpha * sum(abs(B[,i])))
    }
    penalty <- sum(penalty_temp)
  }
  
  return(LL - penalty)
  
}

#' Wrapper function to calculate all log-likelihood function for multivariate elastic net regression

.LL_all_mvr <- function(X, Z, r_i, Y = NULL,  V, W, lambda1, lambda2, alpha1, alpha2,
                       multiply = T, mode = "ps", ratio_max = 20){
  
  
  
  n <- nrow(X)
  p <- ncol(X)
  
  V <- as.matrix(V)
  
  if(is.null(Y)){
    
    mode = "ps"
    
    Y <- X
  }else{
    mode = "reg"
    
    Y <- as.matrix(Y)
  }
  
  nr <- ncol(Y)
  K <- ncol(V)
  q <- ncol(Z)
  
  if(is.null(r_i)){
    
    r_i <- matrix(1/K, nrow = n, ncol = K)
    
  }
  
  if(multiply){
    lambda1 <- lambda1 * n
    lambda2 <- lambda2 * n
  }
  
  if(is.null(dim(lambda1)[2])){

    lambda1 <- matrix(rep(lambda1, nr) ,nrow = K, ncol = nr)
    
  }
  
  sgm <- .est_sgm_K(X, Y, W, r_i)
  
  Z_1 <- cbind(rep(1,n), Z)
  X_1 <- cbind(rep(1,n), X)
  
  pi_est <- prob_multi(Z_1, V)
  
  LL_n <- matrix(0, nrow = n, ncol = K)
  
  penalty_W <- matrix(0, nrow = K, ncol = nr)
  penalty_V <- c()
  
  for(i in 1:K){
    
    LL_n[,i] <- .LLikeli_mvr(X = X, B = W[,,i], Y = Y, sgmY =  sgm[i,], reduce = F, intercept = T)
    for(j in 1:nr){
      penalty_W[i,j] <- (lambda1[i,j] / sgm[i, j]) * (0.5 * (1 - alpha1[i]) * sum(W[2:(p + 1),j,i]^2) + alpha1[i] * sum(abs(W[2:(p + 1),j,i])))
    }
    
  }
  
  penalty_V <- lambda2 * (0.5 * (1 - alpha2) * sum(V[2:(q + 1),]^2) + alpha2 * sum(abs(V[2:(q + 1),])))
  
  LL_comp <- sum(r_i * LL_n)
  PLL_comp <- LL_comp - penalty_V - sum(penalty_W)
  
  LL_diff <- LL_n - LL_n[,1]
  pi_diff <- log(pi_est) - log(pi_est[,1])
  
  total_diff <- LL_diff + pi_diff
  total_diff[total_diff > ratio_max] = ratio_max
  total_diff[total_diff < -ratio_max] = -ratio_max
  
  
  diff_ratio <- log(rowSums(exp(total_diff)))
  
  LL_obs <- sum(LL_n[,1] + log(pi_est[,1]) + diff_ratio)
  PLL_obs <- LL_obs - penalty_V - sum(penalty_W)
  
  return(c(LL_obs, PLL_obs, LL_comp, PLL_comp))
}


.est_sgm <- function(X, Y, B, r_i = NULL){
  
  n <- nrow(X)
  p <- ncol(X)
  
  nr <- ncol(Y)
  
  X_1 <- cbind(rep(1,n), X)
  
  if(is.null(r_i)){
    r_i <- rep(1,n)
  }
  
  sgm <- sqrt(colSums(r_i * (Y - X_1 %*% B)^2)/sum(r_i))
  
  return(sgm)
}

.est_sgm_K <- function(X, Y, B, r_i = NULL){
  
  
  n <- nrow(X)
  p <- ncol(X)
  
  nr <- ncol(Y)
  K <- dim(B)[3]
  
  if(is.null(r_i)){
    r_i <- matrix(1/K, nrow = n, ncol = K)
  }
  
  sgm <- matrix(0, nrow = K, ncol = nr)
  
  for(i in 1:K){
    
    if(sum(r_i[,i]))
    
    sgm[i,] <- sqrt(.est_sgm(X, Y, B = B[,,i], r_i = r_i[,i]))
    
  }
  
  return(sgm)
}

#' objective function for elastic net regression

.obj_elnet <- function(y, X, R = NULL, beta, lambda, alpha, intercept = F, multiply = F){
  
  n <- nrow(X)
  p <- ncol(X)
  
  if(multiply){
    lambda <- lambda * n
  }
  
  if(intercept){
    X <- cbind(rep(1,n),X)
    penalty <-  lambda * (0.5*(1 - alpha)*sum(beta[2:nrow(beta)]^2) + alpha*sum(abs(beta[2:nrow(beta)])))
  }else{
    penalty <-  lambda * (0.5*(1 - alpha)*sum(beta^2) + alpha*sum(abs(beta)))
  }
  
  if(!length(R)){
    R = rep(1, n)
  }else{
    R = as.vector(R)
  }
  
  obj <- 0.5*sum(diag(R) %*% (y - X %*% beta)^2) + penalty
  #  print(obj)
  
  return(obj)
  
}

#' objective function for logistic regression

.obj_logistic <- function(X, beta, y, R = NULL, lambda, alpha,
                         intercept = F, multiply = F, residual_k = 1){
  
  n <- nrow(X)
  p <- ncol(X)
  
  if(multiply){
    lambda <- lambda * n
  }
  
  beta <- as.matrix(beta)
  if(ncol(beta) > 1){
    beta <- as.matrix(beta[,2])
  }
  
  if(!length(R)){
    R = rep(1, n)
  }else{
    R = as.vector(R)
  }
  
  if(intercept){
    X <- cbind(rep(1,n),X)
    penalty <-  lambda * (0.5*(1 - alpha)*sum(beta[2:nrow(beta), ]^2) + alpha*sum(abs(beta[2:nrow(beta), ])))
  }else{
    penalty <-  lambda * (0.5*(1 - alpha)*sum(beta^2) + alpha*sum(abs(beta)))
  }
  
  Likeli <- sum(R*log(.Likeli_binomial(X, beta, y, log = F, residual_k = residual_k)))
  
  binom = - Likeli + penalty
  
  return(binom)
}

#' objective function for sparse multinomial regression

.obj_multi <- function(X, beta, y, R = NULL, lambda, alpha, intercept = F, multiply = F){
  
  n <- nrow(X)
  p <- ncol(X)
  
  if(!length(R)){
    R = rep(1, n)
  }else{
    R = as.vector(R)
  }
  
  y <- as.matrix(y)
  if(ncol(y) == 1){
    y <- cbind( 1- y, y)
  }
  
  if(multiply){
    lambda <- lambda * n
  }
  
  if(intercept){
    X <- cbind(rep(1,n),X)
    penalty <- lambda * (0.5*(1 - alpha)*sum(beta[2:nrow(beta),]^2) + alpha*sum(abs(beta[2:nrow(beta),])))
  }else{
    penalty <- lambda * (0.5*(1 - alpha)*sum(beta^2) + alpha*sum(abs(beta)))
  }
  
  Likeli <- sum(R* log(.Likeli_multi(X, beta, y, log = F)))
  
  multinom = - Likeli + penalty
  
  return(multinom)
}

#' Calculate likelihood for binomial distribution

.Likeli_binomial0 <- function(y, prob, thresh = 1e-9, log = FALSE){
  
  Likelihood = prob^y*(1 - prob)^(1 - y)
  
  if(log){
    
    Likelihood[Likelihood < thresh] = thresh
    
    Likeli = sum(log(Likelihood))
  }else{
    
    Likeli <- Likelihood
    
  }
  
  return(Likeli)
  
}

# Calculate likelihood for binomial distribution with logistic regression input

.Likeli_binomial <- function(X, B, y, thresh = 1e-9, log = FALSE, residual_k = 1){
  
  odd <- X %*% B
  
  prob <- 1/(1 + residual_k * exp(-odd))
  
  Likeli <- .Likeli_binomial0(y, prob, thresh, log = log)
  
  return(Likeli)
  
}

#' Calculate probability for multinomial regression
#' @param X the input data matrix (Need to add intercept term).
#' @param B coefficient of multinomial regression.
#' return the probability for multinomial regression.
#' @export

prob_multi <- function(X, B){
  
  odd <- X %*% B
  
  prob <- t(apply(odd, 1, softmax))
  
  return(prob)
  
}

.Likeli_multi0 <- function(y, prob, thresh = 1e-9, log = FALSE){
  
  Likeli_temp <- prob ^ y
  
  Likelihood = apply(Likeli_temp, 1, prod)
  
  if(log){
    
    Likelihood[Likelihood < thresh] = thresh
    
    Likelihood[is.na(Likelihood)] = 0
    
    Likeli = sum(log(Likelihood))
  }else{
    
    Likeli <- Likelihood
    
  }
  
  return(Likeli)
  
}

#' Calculate likelihood function of multinomial distribution with regression input.

.Likeli_multi <- function(X, B, y, thresh = 1e-9, log = FALSE){
  
  prob <- prob_multi(X, B)
  
  Likeli <- .Likeli_multi0(y, prob, thresh, log= log)
  
  return(Likeli)
}

#' Calculate complete log likelihood function, i.e. the corresponding Q function.

.LogLikeli_complete <- function(X, Z, y, beta_V, beta_W, reduce = FALSE){
  
  n <- nrow(X)
  
  K <- ncol(beta_V)
  
  latent_prob <- prob_multi(Z, beta_V)
  
  L_z <- matrix(0, nrow = n, ncol = K)
  
  for(k in 1:K){
    L_z[,k] <- Likeli_binomial(X = X, y = y, B = beta_W[[k]])
  }
  
  LL_z <- latent_prob*L_z/rowSums(latent_prob*L_z)
  
  temp <- latent_prob*L_z
  
  temp[temp < 1e-9] <- 1e-9
  
  LL_f <- log(temp)
  
  LLikeli <- rowSums(LL_z*LL_f)
  
  if(reduce){
    LLikeli <- sum(LLikeli)
  }
  
  return(LLikeli)
}

#' Calculate the observed likelihood function, i.e. likelihood function for mixture distribution.

.Likeli_obs <- function(X, Z, y, beta_V, beta_W, log = FALSE, thresh = 1e-6){
  
  n <- nrow(X)
  
  K <- ncol(beta_V)
  
  latent_prob <- prob_multi(Z, beta_V)
  
  sub_prob <- matrix(0, nrow = n, ncol = K)
  
  for(i in 1:K){
    
    sub_prob[,i] <- 1/(1 + exp(-(X %*% beta_W[[i]])))
    
  }
  
  prob <- rowSums(latent_prob * sub_prob)
  
  Likelihood <- prob^y*(1 - prob)^(1 - y)
  
  if(log){
    
    Likelihood[Likelihood < thresh] = thresh
    
    Likeli = sum(log(Likelihood))
  }else{
    
    Likeli <- Likelihood
    
  }
  
  return(Likeli)
}

#' Calculate penalized likelihood function.

.PLL <- function(X, Z, y, beta_V, beta_W, version = "complete",
                lambda1 = NULL, lambda2 = NULL, alpha1 = 1, alpha2 = 1, multiply = F){
  
  n <- nrow(X)
  
  K <- ncol(beta_V)
  
  p <- length(beta_W[[1]])
  q <- nrow(beta_V)
  
  if(multiply){
    lambda1 <- lambda1 * n
    lambda2 <- lambda2 * n
  }
  
  penalty1 <- 0
  
  if(length(lambda1) == 1){
    lambda1 <- rep(lambda1, K)
  }
  
  for(i in 1:K){
    penalty1_temp <- lambda1[i] * (0.5 * (1 - alpha1[i]) * sum(beta_W[[i]][2:p]^2) + alpha1[i] * sum(abs(beta_W[[i]][2:p])))
    penalty1 <- penalty1 + penalty1_temp
  }
  
  penalty2 <- lambda2 * (0.5 * (1 - alpha2) * sum(beta_V[2:q,]^2) + alpha2 * sum(abs(beta_V[2:q,])))
  
  if(version == "obs"){
    Likeli <- Likeli_obs(X, Z, y, beta_V, beta_W)
    
    PLLikeli <- sum(log(Likeli)) - penalty1 - penalty2
    
  }else{
    LLikeli <- LogLikeli_complete(X, Z, y, beta_V, beta_W)
    
    PLLikeli <- sum(LLikeli) - penalty1 - penalty2
  }
  
  
  return(PLLikeli)
}

#' Q function for given current parameters.

.Q_func <- function(X, Z, y, beta_V, beta_W, beta_V0, beta_W0,
                   lambda1 = NULL, lambda2 = NULL, alpha1 = 1, alpha2 = 1,
                   reduced = FALSE, multiply = F){
  
  n <- nrow(X)
  
  K <- ncol(beta_V)
  p <- ncol(X)
  q <- nrow(beta_V)
  
  penalty1 <- 0
  penalty10 <- 0
  
  if(length(lambda1) == 1){
    lambda1 <- rep(lambda1, K)
  }
  
  if(multiply){
    lambda1 <- lambda1 * n
    lambda2 <- lambda2 * n
  }
  
  for(i in 1:K){
    penalty1_temp <- lambda1[i] * (0.5 * (1 - alpha1[i]) * sum(beta_W[[i]][2:p]^2) + alpha1[i] * sum(abs(beta_W[[i]][2:p])))
    penalty1 <- penalty1 + penalty1_temp
    penalty10_temp <- lambda1[i] * (0.5 * (1 - alpha1[i]) * sum(beta_W0[[i]][2:p]^2) + alpha1[i] * sum(abs(beta_W0[[i]][2:p])))
    penalty10 <- penalty10 + penalty10_temp
  }
  
  penalty2 <- lambda2 * (0.5 * (1 - alpha2) * sum(beta_V[2:q,]^2) + alpha2 * sum(abs(beta_V[2:q,])))
  penalty20 <- lambda2 * (0.5 * (1 - alpha2) * sum(beta_V0[2:q,]^2) + alpha2 * sum(abs(beta_V0[2:q,])))
  
  latent_prob <- prob_multi(Z, beta_V0)
  pi <- prob_multi(Z, beta_V)
  
  L_z0 <- matrix(0, nrow = n, ncol = K)
  L_z <- matrix(0, nrow = n, ncol = K)
  
  for(k in 1:K){
    L_z0[,k] <- Likeli_binomial(X = X, y = y, B = beta_W0[[k]])
    L_z[,k] <- Likeli_binomial(X = X, y = y, B = beta_W[[k]])
  }
  
  LL_z <- latent_prob*L_z0/rowSums(latent_prob*L_z0)
  LL_Ent <- sum(LL_z * log(pi))
  LL_Ent0 <- sum(LL_z * log(latent_prob))
  
  LL_f <- log(pi*L_z)
  LL_f0 <- log(latent_prob *L_z0)
  
  LLikeli <- colSums(LL_z*LL_f)
  LLikeli0 <- colSums(LL_z*LL_f0)
  
  Q = sum(LLikeli)
  Q0 = sum(LLikeli0)
  
  if(reduced){
    return(Q- penalty1 - penalty2)
  }else{
    return(list(LL_sub = LLikeli, LL_sub0 = LLikeli0,
                LL_Ent = LL_Ent, LL_Ent0 = LL_Ent0,
                P1 = penalty1, P2 = penalty2,
                Q = Q, PQ = Q - penalty1 - penalty2,
                Q0 = Q0, PQ0 = Q0 - penalty10 - penalty20))
  }
  
}

#' Function to calculate loglikelihood function including observed loglikelihood,
#' penalized loglikehood, complete loglikelihood and penalized complete loglikelihood.

.LL_all <- function(X, Z, y, beta_V, beta_W, lambda1, lambda2, alpha1, alpha2, multiply = F){
  
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  
  y <- as.matrix(y)
  
  X_1 <- cbind(rep(1,nrow(X)), X)
  Z_1 <- cbind(rep(1,nrow(Z)), Z)
  
  LL_obs <- .Likeli_obs(X = X_1, Z = Z_1, y = y,
                       beta_V = beta_V, beta_W = beta_W,
                       log = TRUE)
  PLL_obs <- .PLL(X = X_1, Z = Z_1, y = y,
                 beta_V = beta_V, beta_W = beta_W, version = "obs",
                 lambda1 = lambda1, lambda2 = lambda2,
                 alpha1 = alpha1, alpha2 = alpha2, multiply = multiply)
  
  LL_complete <- .LogLikeli_complete(X = X_1, Z = Z_1, y = y,
                                    beta_V = beta_V, beta_W = beta_W,
                                    reduce = TRUE)
  PLL_complete <- .PLL(X = X_1, Z = Z_1, y = y,
                      beta_V = beta_V, beta_W = beta_W,version = "complete",
                      lambda1 = lambda1, lambda2 = lambda2,
                      alpha1 = alpha1, alpha2 = alpha2, multiply = multiply)
  
  
  
  LL_res <- c(LL_obs, PLL_obs, LL_complete, PLL_complete)
  
  names(LL_res) <- c("LL_obs", "PLL_obs", "LL_complete", "PLL_complete")
  
  return(LL_res)
}

#' Function to calculate loglikelihood function with list input.

.LL_all_mat <- function(X_list, Z, y, beta_V, beta_W, lambda1_mat, lambda2, alpha1_mat, alpha2, multiply = F){
  
  l <- length(X_list)
  
  LL_res <- matrix(0, nrow = l, ncol = 4)
  rownames(LL_res) <- names(X_list)
  
  for(i in 1:l){
    
    X_temp <- X_list[[i]]
    y_temp <- y[,i]
    
    W_temp <- mat2list(beta_W[[i]])
    
    LL_res[i,] <- LL_all(X_temp, Z, y_temp, beta_V, W_temp, lambda1_mat[i,], lambda2, alpha1_mat[i,], alpha2, multiply = multiply)
    
  }
  
  colnames(LL_res) <- c("LL_obs", "PLL_obs", "LL_complete", "PLL_complete")
  
  return(LL_res)
  
}

#' Calculate the probability for each latent class with multiple input.

.r_i_arr <- function(X_list, Z, y, pi_old, beta_W, beta_V, r_i_eps = 1e-6){
  
  
  l <- length(X_list)
  n <- nrow(Z)
  K <- ncol(beta_V)
  
  r_i <- array(0, dim = c(n, K, l))
  
  for(i in 1:l){
    
    X_temp <- X_list[[i]]
    X_1 <- cbind(rep(1,n), X_temp)
    y_1 <- y[,i]
    
    Likeli_W <- purrr::partial(Likeli_binomial, X = X_1, y = y_1, log = FALSE)
    W_old <- mat2list(beta_W[[i]])
    Likeli_old <- sapply(W_old, Likeli_W)
    
    r_i_new0 <- pi_old * Likeli_old
    r_i_new_s <- rowSums(r_i_new0)
    r_i_new <- r_i_new0
    
    r_i_new[r_i_new_s > r_i_eps, ] <- r_i_new0[r_i_new_s > r_i_eps, ] / r_i_new_s[r_i_new_s > r_i_eps]
    r_i_new[r_i_new_s <= r_i_eps, ] <- 1/K
    
    r_i[,,i] <- r_i_new
  }
  
  return(r_i)
}


