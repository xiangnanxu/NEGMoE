# Inner loop for Elastic net regression updating.

.Elnet_step <- function(X, y, R = NULL, beta_old, lambda, alpha = 1, intercept = TRUE){

  n <- nrow(X)
  p <- ncol(X)

  if(!length(R)){
    R = rep(1, n)
  }else{
    R = as.vector(R)
  }

  for(i in 1:p){

    Xi <- X[,i]
    X_i <- X[,-i]
#    betai <- beta_old[i,]
#    beta_i <- beta_old[-i,]
    beta_i <- beta_old
    beta_i[i,] = 0
    Xi_normsq <- sum(Xi^2)
    Xi_normsq_w <- sum(R * Xi^2)

    X_i <- as.matrix(X_i)

    if(intercept){
      if(i == 1){

#        beta_temp <- sum (R*(y - X_i %*% beta_i))/ sum(R)
        beta_temp <- sum (R*(y - X %*% beta_i))/ sum(R)
        beta_temp_thresh <- beta_temp

      }else{
#        beta_temp <- sum(R * Xi * (y - X_i %*% beta_i))
        beta_temp <- sum(R * Xi * (y - X %*% beta_i))
        beta_temp_thresh <- soft_thresh(beta_temp, lambda * alpha)
        beta_temp_thresh <- beta_temp_thresh / (Xi_normsq_w + lambda * (1 - alpha))
      }
    }else{
      beta_temp <- sum(R * Xi * (y - X_i %*% beta_i))
      beta_temp_thresh <- soft_thresh(beta_temp, lambda * alpha)
      beta_temp_thresh <- beta_temp_thresh / (Xi_normsq_w + lambda * (1 - alpha))
    }
    beta_old[i,] <- beta_temp_thresh

  }

  return(beta_old)
}

# Coordinate descent method solving elastic net regression.

.Elnet <- function(X, y, lambda, alpha = 1, R = NULL, beta_init = NULL,
                   itmax = 1e3, eps = 1e-9, eps_var = 1e-5, intercept = TRUE){

  n <- nrow(X)
  p <- ncol(X)

  y[is.infinite(y)] <- 0
  y[is.na(y)] <- 0

  if(!length(R)){
    R = rep(1, n)
  }else{
    R = as.vector(R)
  }

  X_1 <- cbind(rep(1, n), X)

  lambda = lambda * n

  if(!length(beta_init)){

    beta <- matrix(0, nrow = (p + 1), ncol = 1)

  }else{

    beta <- beta_init

  }

  beta_old <- beta

  if(var(y) < eps_var){

    rownames(beta) <- c("intercept", paste("V",1:p,sep = ""))

    return(beta)

  }else{

    it = 0

    repeat{
      beta <- .Elnet_step(X = X_1, y = y, R = R, beta_old = beta_old,
                         lambda = lambda, alpha = alpha, intercept = intercept)

      cond1 <- (sum(abs(beta - beta_old)) < eps)
      cond2 <- (it > itmax)

      if(cond1|cond2){

        break()

      }

      it = it + 1

      beta_old <- beta

    }

    rownames(beta) <- c("intercept", paste("V", 1:p,sep = ""))
  }
  return(beta)
}

#' Wrapper function for elastic net regression with glmnet.

.Elnet_glmnet <- function(X, y, lambda, alpha = 1, R = NULL, beta_init = NULL,
                         itmax = 1e3, eps = 1e-3, eps_var = 1e-5, intercept = TRUE){

  n <- nrow(X)
  p <- ncol(X)

  y[is.infinite(y)] <- 0
  y[is.na(y)] <- 0

  lambda <- lambda

  if(!length(R)){
    R = rep(1, n)
  }else{
    R = as.vector(R)
  }

  if(var(y) < eps_var){

    beta <- matrix(0, nrow = (p + 1), ncol = 1)

    rownames(beta) <- c("intercept", paste("V",1:p,sep = ""))

    return(beta)

  }else{

    beta <- glmnet::glmnet(x = X, y = y, family = "gaussian",
                           weights = as.vector(R), lambda = lambda * (n / sum(R)),
                           alpha = alpha, standardize = F)

    return(as.matrix(coef(beta)))
  }


}
