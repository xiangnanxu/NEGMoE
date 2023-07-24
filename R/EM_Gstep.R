# Update function for CEM algorithm
.CEM_Gstep <- function(X, Z, r_i_new, lambda1, lambda2, alpha2, beta_max = 10, itmax0 = 1e2,
                     V_old = NULL, Theta_old = NULL, mu_old = NULL,
                     method = "icov", thresh_cor = 0.3, thresh_eig = 0.1){
  
  n <- nrow(X)
  p <- ncol(X)
  
  q <- ncol(Z)
  K <- ncol(r_i_new)
  
  X_1 <- cbind(rep(1,nrow(X)), X)
  Z_1 <- cbind(rep(1,nrow(Z)), Z)
  
  r_i <- hard_prob(r_i_new)
  

  Theta_new <- array(0, dim = c(p,p,K))
  mu_new <- matrix(0, nrow = K, ncol = p)
  S_new <- array(0, dim = c(p,p,K))
  
  for(i in 1:K){
    
    X_k <- X[r_i[,i] == 1,]
    
    mu_k <- colMeans(X_k)
    mu_new[i,] <- mu_k

    if(method == "icov"){
      
      Theta_k <- GLasso(X = X_k, mu = mu_k, lambda = lambda1[i], Theta_init = Theta_old[,,i], itmax = itmax0)$Theta
      
      S_k <- chol_solve(Theta_k)
      
    }else{
      
      S_k <- cov(X_k)
      S_k <- proj_sd_thresh(S_k, thresh_cor = thresh_cor, thresh_eig)
      
      Theta_k <- chol_solve(S_k)
      
    }
    
    
    Theta_new[,,i] <- Theta_k
    S_new[,,i] <- S_k
  }
  
  V_new <- try(sMulti_glmnet(X = Z, y = r_i, lambda = lambda2, alpha = alpha2, beta_max = beta_max))
  if("try-error" %in% class(V_new)){
    if(length(V_old)){
      V_new <- .sMulti(X = Z, y = r_i, lambda = lambda2, alpha = alpha2, beta_max = beta_max, beta_init = V_old)
    }else{
      V_new <- .sMulti(X = Z, y = r_i, lambda = lambda2, alpha = alpha2, beta_max = beta_max)
    }
    
  }
  
  return(list(Theta_new = Theta_new, S_new = S_new, mu_new = mu_new, V_new = V_new))
  
}

# Update function for SEM algorithm
.SEM_Gstep <- function(X, Z, r_i_new, lambda1, lambda2, alpha2, beta_max = 10, itmax0 = 1e2,
                     V_old = NULL, Theta_old = NULL, mu_old = NULL,
                     method = "icov", thresh_cor = 0.3, thresh_eig = 0.1){
  
  K <- ncol(r_i_new)
  
  it <- 0
  repeat{
    
    r_i <- sample_prob(r_i_new)
    
    if(class_check(r_i)){
      break()
    }
    
    it <- it + 1
    if(it > 100){
      a <- rep(1, K)
      r_i_0 <- rdirichlet(n, a)
      r_i <- sample_prob(r_i_0)
      break()
    }
  }
  
  return(.CEM_Gstep(X = X, Z = Z, r_i_new = r_i, lambda1 = lambda1, itmax0 = itmax0,
                  lambda2 = lambda2, alpha2 = alpha2, beta_max = beta_max,
                  V_old = V_old, Theta_old = Theta_old, mu_old = mu_old,
                  method = method, thresh_cor = thresh_cor, thresh_eig = thresh_eig))
}

# Update function for SAEM algorithm
.SAEM_Gstep <- function(X, Z, r_i_new, lambda1, lambda2, itmax0 = 1e2,
                      alpha2, itmax = 1e3, beta_max = 10,
                      V_old, Theta_old, mu_old,
                      method = "icov", thresh_cor = 0.3, thresh_eig = 0.1){
  
  K <- ncol(r_i_new)
  n <- nrow(X)
  
  it <- 0
  repeat{
    
    r_i <- sample_prob(r_i_new)
    
    if(class_check(r_i)){
      break()
    }
    
    it <- it + 1
    if(it > 100){
      a <- rep(1, K)
      r_i_0 <- rdirichlet(n, a)
      r_i <- sample_prob(r_i_0)
      break()
    }
  }
  
  r_i_new <- (1 - gamma) * r_i_new + gamma * r_i
  
  return(.BEM_Gstep(X = X, Z = Z, r_i_new = r_i_new,
                  lambda1, lambda2, alpha2, itmax0 = itmax0,
                  V_old, Theta_old, mu_old,
                  method = method, thresh_cor = thresh_cor, thresh_eig = thresh_eig))
}

# Update function for EM algorithm
.BEM_Gstep <- function(X, Z, r_i_new, lambda1, lambda2,
                     alpha2, itmax0 = 10, beta_max = 10,
                     V_old, Theta_old, mu_old,
                     method = "icov", thresh_cor = 0.3, thresh_eig = 0.1){
  
  
  n <- nrow(X)
  p <- ncol(X)
  
  q <- ncol(Z)
  K <- ncol(r_i_new)
  
  X_1 <- cbind(rep(1,nrow(X)), X)
  Z_1 <- cbind(rep(1,nrow(Z)), Z)
  
  Theta_new  = array(0, dim = c(p,p,K))
  S_new = array(0, dim = c(p,p,K))
  mu_new <- matrix(0, nrow = K, ncol = p)
  
  for(i in 1:K){
    
    mu_new[i,] <- colSums(r_i_new[,i] * X) / sum(r_i_new[,i])

    if(method == "icov"){
      
      Theta_k <- GLasso(X, mu_new[i,], lambda = lambda1[i], weight = r_i_new[,i], itmax = itmax0, Theta_init = Theta_old[,,i])$Theta
      
      S_k <- chol_solve(Theta_k)
      
    }else{
      
      S_k <- (t(X - mu_new[i,]) %*% diag(r_i_new[,i]) %*% (X - mu_new[i,]))/sum(r_i_new[,i])
      S_k <- proj_sd_thresh(S_k, thresh_cor = thresh_cor, thresh_eig = thresh_eig)
      
      Theta_k <- chol_solve(S_k)
    }
    
    Theta_new[,,i] <- Theta_k
    S_new[,,i] <- S_k
    
  }
  
  V_new <- try(.sMulti(X = Z, y = r_i_new, lambda = lambda2, alpha = alpha2, beta_max = beta_max, version = "glmnet"))
  if("try-error" %in% class(V_new)){
    if(length(V_old)){
      V_new <- .sMulti(X = Z, y = r_i_new, lambda = lambda2, alpha = alpha2, beta_max = beta_max, beta_init = V_old, version = "my")
    }else{
      V_new <- .sMulti(X = Z, y = r_i_new, lambda = lambda2, alpha = alpha2, beta_max = beta_max, version = "my")
    }
    
  }
  
  return(list(Theta_new = Theta_new, S_new =S_new, mu_new = mu_new, V_new = V_new))
}

# Wrapper function for EM step update
.EM_Gstep <- function(X, Z, lambda1, lambda2, alpha2, gamma = 1, beta_max = 10, ratio_max = 20,
                    V_old, mu_old, Theta_old, itmax0 = 1e2, backtracking = "comp", EM_opt = "EM",
                    method = "icov", thresh_cor = 0.3, thresh_eig = 0.1, adapt = T){
  n <- nrow(X)
  p <- ncol(X)
  
  q <- ncol(Z)
  K <- ncol(V_old)
  
  X_1 <- cbind(rep(1,nrow(X)), X)
  Z_1 <- cbind(rep(1,nrow(Z)), Z)
  
  
  pi_old <- prob_multi(Z_1, V_old)
  
  Likeli_G <- matrix(0, nrow = n, ncol = K)
  for(i in 1:K){
    Likeli_G[,i] <- .LLglasso(X = X, mu = mu_old[i,], Theta = Theta_old[,,i], reduce = F)
    
  }
  
  pi_ratio <- log(pi_old) - log(pi_old[,1])
  Likeli_ratio <- Likeli_G - Likeli_G[,1]
  
  total_ratio <- pi_ratio + Likeli_ratio
  
  total_ratio[total_ratio > ratio_max] = ratio_max
  total_ratio[total_ratio < -ratio_max] = -ratio_max
  
  r_i_new <- t(apply(total_ratio, 1, softmax))
  
  if(adapt){
    pen_fac <- rep(1, K)
  }else{
    
    psu_sz <- colSums(r_i_new)
    psu_sz[psu_sz < 1] = 1
    
    pen_fac <- n/(psu_sz)
  }
  
  
  lambda1_fac <- lambda1 * pen_fac
  
  if(EM_opt %in% c("SEM", "SAEM", "EM")){
    V_old <- clipping(V_old, beta_max, -beta_max)
  }
  
  
  if(EM_opt == "CEM"){
    step_func <- purrr::partial(.CEM_Gstep, Theta_old = Theta_old, mu_old = mu_old, V_old = V_old,
                                method = method, thresh_cor = thresh_cor, thresh_eig = thresh_eig)
  }else if(EM_opt == "SEM"){
    step_func <- purrr::partial(.SEM_Gstep, Theta_old = Theta_old, mu_old = mu_old, V_old = V_old,
                                method = method, thresh_cor = thresh_cor, thresh_eig = thresh_eig)
  }else if(EM_opt == "EM"){
    step_func <- purrr::partial(.BEM_Gstep, itmax = itmax0, beta_max = beta_max,
                                Theta_old = Theta_old, mu_old = mu_old, V_old = V_old,
                                method = method, thresh_cor = thresh_cor, thresh_eig = thresh_eig)
  }else if(EM_opt == "SAEM"){
    step_func <- purrr::partial(.SAEM_Gstep, gamma = gamma, itmax = itmax0, beta_max = beta_max,
                                Theta_old = Theta_old, mu_old = mu_old, V_old = V_old,
                                method = method, thresh_cor = thresh_cor, thresh_eig = thresh_eig)
  }
  
  EM_step_res <- step_func(X = X, Z = Z, r_i_new = r_i_new,
                           lambda1 = lambda1_fac, lambda2 = lambda2, alpha2 = alpha2)
  
  
  Theta_new <- EM_step_res$Theta_new
  S_new <- EM_step_res$S_new
  V_new <- EM_step_res$V_new
  mu_new <- EM_step_res$mu_new
  

    
  if(backtracking == "none"){
    
    Theta_new <- Theta_new
    S_new <- S_new
    V_new <- V_new
    
  }else{
    
    obj_mult <- purrr::partial(.obj_multi, X = Z, y = r_i_new, 
                               lambda = lambda2, alpha = alpha2, intercept = T, multiply = T)
    
    V_btr <- .backtracking_f(obj_mult, x_old = V_old, x_new = V_new, mm = "min")
    Theta_btr <- array(0, dim = c(p, p, K))
    S_btr <- array(0, dim = c(p, p, K))
    
    for(i in 1:K){
      
      if(method == "icov"){
        
        obj_Glasso <- purrr::partial(.PLLglasso, X = X, mu = mu_new[i,], weight = r_i_new[,i], lambda = lambda1_fac[i])
        
        Theta_temp <- .backtracking_f(obj_Glasso, x_old = Theta_old[,,i], x_new = Theta_new[,,i], mm = "max")
        
        Theta_temp <- (Theta_temp + t(Theta_temp))/2
        Theta_temp <- proj_sd(Theta_temp, thresh = thresh_eig)
        
        S_temp <- chol_solve(Theta_temp)
        
      }else{
        
        S_old <- chol_solve(Theta_old[,,i])
        
        obj_Glasso <- purrr::partial(.LLglasso, X = X, mu = mu_new[i,], weight = r_i_new[,i], type = "cov")
        
        S_temp <- .backtracking_f_proj(obj_Glasso, x_old = S_old, x_new = S_new[,,i], mm = "max", 
                                      thresh_cor = thresh_cor, thresh_eig = thresh_eig)
        
        
        S_temp <- (S_temp + t(S_temp))/2
        S_temp <- proj_sd_thresh(S_temp, thresh_cor = thresh_cor, thresh_eig = thresh_eig)
        
        Theta_temp <- chol_solve(S_temp)
      }
      
      
      
      Theta_btr[,,i] <- Theta_temp
      
      S_btr[,,i] <- S_temp
    }
    
    V_new <- hard_thresh(V_btr)
    Theta_new <- Theta_btr
    S_new <- S_btr
    
  }
  
  LL_new <- LL_Glasso(X = X, Z = Z, beta_V = V_new, Theta = Theta_new, mu = mu_new, 
                      r_i = r_i_new, lambda1 = lambda1, lambda2 = lambda2)

  
  return(list(Theta_new = Theta_new, S_new = S_new, mu_new = mu_new,
              V_new = V_new, r_i = r_i_new, pi = pi_old, LL_new = LL_new))
}
