###################################################################
# function to fitting NEGMoE using Graphical Lasso
###################################################################

.EM_GMoE0 <- function(X, Z, Y = NULL, K = 2, lambda1 = 0.2, lambda2 = 0.03,alpha1 = 1,
                      alpha2 = 1, stop = "all", gamma_func = inv_func, beta_max = 10,
                      solver = "glmnet", stop_eps1 = 1e-2, stop_eps2 = 1e-2,
                      itmax = 1e2, itmin = 3, itmax0 = 10, itmax1 = 1e2, itmax2 = 1e2,
                      W_init = NULL, V_init = NULL, EM_opt = "EM", backtracking = "comp",
                      adapt = T, init = "kmeans",ratio_max = 20, method = "icov",
                      thresh_cor = 0.3, thresh_eig = 0.1, verbose = TRUE){
  
  if(!(EM_opt %in% c("CEM", "GEM", "SEM", "SAEM", "EM"))){
    message("Algorithm not implemented, EM as default")
    EM_opt = "EM"
  }
  
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  
  V_name <- colnames(Z)
  Theta_name <- colnames(X)
  
  if(is.null(V_name)){
    V_name <- c("intercept", paste("V",1:q, sep = ""))
  }else{
    V_name <- c("intercept", V_name)
  }
  
  if(is.null(Theta_name)){
    Theta_name <- c(paste("V",1:p, sep = ""))
  }else{
    Theta_name <- c(Theta_name)
  }
  

  X_1 <- cbind(rep(1,nrow(X)), X)
  Z_1 <- cbind(rep(1,nrow(Z)), Z)
  
  if(!(length(V_init))){
    if(init == "rand"){
      a <- rep(1, K)
      r_i <- rdirichlet(n, a)
      r_i <- hard_prob(r_i)
      
    }else{
      temp <- kmeans(Z, centers = K)$cluster
      r_i <- onehot(temp)
    }
    
  }else{
    r_i <- prob_multi(Z_1, V_init)
    r_i_0 <- hard_prob(r_i)
    
    if(!class_check(r_i_0)){
      r_i <- sample_prob(r_i)
    }else{
      r_i <- r_i_0
    }
  }
  
  if(!(length(V_init))){
    V_old <- sMulti_glmnet(Z, r_i, lambda = lambda2, alpha = alpha2, beta_max = beta_max)
  }else{
    V_old <- V_init
  }
  
  lambda1 <- ext_param(lambda1, K)
  
  S_0 <- cov(X)
  S_old <- array(0, dim = c(p,p,K))
  
  if(is.null(W_init)){
    Theta_0 <- array(0, dim = c(p,p,K))
    Theta_old <- array(0, dim = c(p,p,K))
    
    mu_old <- matrix(0, nrow = K, ncol = p)
    
    for(i in 1:K){
      
      mu_old[i,] <- colSums(r_i[,i] * X) / sum(r_i[,i])
      S_temp <- (t(X - mu_old[i,]) %*% diag(r_i[,i]) %*% (X - mu_old[i,]))/sum(r_i[,i])
      
      if(method == "icov"){
        
        if(all(lambda1 == 0)){
          Theta_temp <- chol_solve(S_old[,,i])
        }else{
          Theta_temp <- GLasso(X = X, mu = mu_old[i,], lambda = lambda1[i],
                               weight = r_i[,i], itmax = itmax0)$Theta
        }
      }else{
        S_temp <- proj_sd_thresh(S_temp, thresh_cor = thresh_cor, thresh_eig = thresh_eig)
        Theta_temp <- chol_solve(S_temp)
      }
      
      S_old[,,i] <- S_temp
      Theta_old[,,i] <- Theta_temp
      
    }
    
    #  print(r_i)
    if(method == "icov"){
      if(all(lambda1 == 0)){
        Theta_0 <- chol_solve(S_0)
      }else{
        Theta_0 <- GLasso(X = X, lambda = mean(lambda1), itmax = itmax0)$Theta
      }
    }else{
      S_0 <- proj_sd_thresh(S_0, thresh_cor = thresh_cor, thresh_eig = thresh_eig)
      Theta_0 <- chol_solve(S_0)
    }
  }else{
    mu_old <- matrix(0, nrow = K, ncol = p)
    
    for(i in 1:K){
      
      mu_old[i,] <- colSums(r_i[,i] * X) / sum(r_i[,i])
      Theta_old = W_init
      Theta_0 = apply(Theta_old, c(1,2), mean)
    }
  }
  
  
  
  #  W_1 <- sMulti(X, y, lambda = mean(lambda1), alpha = alpha1[1], version = "my")
  
  LL0 <- .LLglasso(X = X, Theta = Theta_0)
  
  PLL0 <- .PLLglasso(X = X, Theta = Theta_0, lambda = mean(lambda1))
  
  #  LL1 <- -obj_logistic(y = y, X = X, beta = W_1,
  #                       lambda = 0, alpha = alpha1[1], intercept = TRUE)
  
  #  PLL1 <- -obj_logistic(y = y, X = X, beta = W_1,
  #                        lambda = mean(lambda1), alpha = alpha1[1], intercept = TRUE, multiply = TRUE)
  
  it = 1
  
  obj <- matrix(0, nrow = (itmax + 1), ncol = 4)
  stop_cond <- matrix(F, nrow = (itmax + 1), ncol = 3)
  obj[1, ] <- c(LL0, PLL0, LL0, PLL0)
  
  repeat{
    
    gamma <- gamma_func(it)
    
    EM_step_res <- .EM_Gstep(X, Z, lambda1, lambda2, alpha2, gamma = gamma,
                             beta_max = beta_max, ratio_max = 20,
                             V_old, mu_old, Theta_old, itmax0 = itmax0,
                             backtracking = backtracking, EM_opt = EM_opt, adapt = adapt,
                             method = method, thresh_cor = thresh_cor, thresh_eig = thresh_eig)
    
    if("try-error" %in% class(EM_step_res)){
      
      message("Early stop! Please check the parameter setting.")
      obj <- obj[1:it,]
      obj <- matrix(obj, ncol = 4)
      colnames(obj) <- c("LL_obs", "PLL_obs", "LL_complete", "PLL_complete")
      stop_cond <- stop_cond[1:it,]
      stop_cond <- matrix(stop_cond, ncol = 3)
      colnames(stop_cond) <- c("it", "cov", "Loss")
      colnames(V_new) <- NULL
      rownames(V_new) <- c("intercept", paste("V",1:q, sep = ""))
      
      return(list(Theta = EM_step_res$Theta_new, V = EM_step_res$V_new, mu = EM_step_res$mu_new,
                  obj = obj, stop_cond = stop_cond, LL0 = LL0, PLL0 = PLL0))
    }
      
      
    V_new <- EM_step_res$V_new
    
    Theta_new <- EM_step_res$Theta_new
    
    S_new <- EM_step_res$S_new
    
    mu_new <- EM_step_res$mu_new
    
    r_i <- EM_step_res$r_i
    pi <- EM_step_res$pi
    
    obj[it + 1,] <- EM_step_res$LL_new
    
    cond_it <- .stop_it(it, itmax = itmax)
    cond_V <- .stop_cov(beta_old = V_old, beta_new = V_new, eps = stop_eps1)
    cond_Theta <- .stop_cov(beta_old = Theta_old, beta_new = Theta_new, eps = stop_eps1)
    cond_S <- .stop_cov(beta_old = S_old, beta_new = S_new, eps = stop_eps1)
    cond_loss <- .stop_loss(obj[,2], it + 1, eps = stop_eps2)
    
    if(method == "icov"){
      stop_cond[it + 1,] <- c(cond_it, cond_V & cond_Theta, cond_loss)
    }else{
      stop_cond[it + 1,] <- c(cond_it, cond_V & cond_S, cond_loss)
    }
    
    
    if(stop == "all"){
      if((it >= itmin)){
        if(all(stop_cond[it + 1,2:3])){
          break()
        }
        else if(cond_it){
          message("Maximium iteration reached!")
          break()
        }
      }
      
    }else{
      if((it >= itmin)){
        if(any(stop_cond[it + 1,2:3])){
          
          break()
        }else if(cond_it){
          message("Maximium iteration reached!")
          break()
        }
      }
    }
    
    V_old <- V_new
    Theta_old <- Theta_new
    S_old <- S_new
    mu_old <- mu_new
    
    if(verbose){
      if(it == 1){
        print(paste0("Fitting NEGMoE (Graphical Lasso)...."))
        print(paste0("it: ",0,", PLL0:",PLL0))
        print(paste0("it: ", it ,", PLL0:",EM_step_res$LL_new[2]))
      }else{
        print(paste0("it: ", it ,", PLL0:",EM_step_res$LL_new[2]))
      }
    }
    
    it = it + 1
  }

  obj <- obj[1:it + 1,]
  obj <- matrix(obj, ncol = 4)
  colnames(obj) <- c("LL_obs", "PLL_obs", "LL_complete", "PLL_complete")
  stop_cond <- stop_cond[1:it + 1,]
  stop_cond <- matrix(stop_cond, ncol = 3)
  colnames(stop_cond) <- c("it", "cov", "Loss")
  colnames(V_new) <- NULL
  rownames(V_new) <- V_name
  rownames(Theta_new) <- Theta_name
  colnames(Theta_new) <- Theta_name
  rownames(mu_new) <- NULL
  colnames(mu_new) <- Theta_name
  
  return(list(Theta = Theta_new, S = S_new, mu = mu_new, V= V_new, r_i = r_i, pi = pi,
              obj = obj, stop_cond = stop_cond, LL0 = LL0, PLL0 = PLL0,
              lambda1 = lambda1, lambda2 = lambda2, alpha2 = alpha2))
}