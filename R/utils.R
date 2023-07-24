
#' Extenion of parameters. To make the number of penalty the same number of components
#' @param lambda Penalty parameters.
#' @param K Number of components.
#' @return A matrix suits for NEGMoE.

ext_param <- function(lambda, K){

  if(length(lambda)  < K){
    lambda <- c(lambda, rep(lambda[length(lambda)], (K - length(lambda))))
  }else if(length(lambda) > K){
    lambda <- lambda[1:K]
  }

  return(lambda)
}

#' Matrix version of Extension of parameters.
#' @param lambda Penalty parameters.
#' @param K Number of components.
#' @param l Number of layers
#' @return A list suits for NEGMoE.
#'
ext_param_mat <- function(lambda, K, l){


  if(length(lambda) == 1){

    lambda1_list <- ext_param(lambda, K * l)

    lambda1_list <- matrix(lambda1_list, nrow = l, ncol = K)

  }else if(length(lambda) == l){

    lambda1_list <- rep(lambda, K)

    lambda1_list <- matrix(lambda1_list, nrow = l, ncol = K)

  }else if(length(lambda) == K){

    lambda1_list <- rep(lambda, each = l)

    lambda1_list <- matrix(lambda1_list, nrow = l, ncol = K)

  }else if(length(lambda) == K*l){
    lambda1_list <- matrix(lambda)

  }else{
    lambda1_list <- rep(lambda1_list, K*l)

    lambda1_list <- lambda1_list[1:K*l]

    lambda1_list <- matrix(lambda1_list, nrow = l, ncol = K)
  }

  return(lambda1_list)
}

#' Generation data follows Dirichlet distribution.
#' @param n Number of samples
#' @param a A vector of length p for Dirichlet distribution parameters.
#' @return A matrix of n x p each row iid follows Dir(a)
#'
rdirichlet <- function(n,a)
{
  if(length(n) > 1 || length(n) < 1 || n < 1) stop("n must be a single positive integer value")
  if(length(a) < 2) stop("a must be a vector of numeric value")
  n <- floor(n)
  l<-length(a);
  x<-matrix(rgamma(l*n,a),ncol=l,byrow=TRUE);
  return(x/rowSums(x) );
}

#' Generation data follows Dirichlet-Multinomial distribution
#' @param n Number of samples
#' @param p Number of variables.
#' @param lib_size library size for each row.(rowSums for each samples)
#' @param pi probability for multinomial distribution. If not provide, will generate by rnorm/sum(rnorm)
#' @param over_disp Overdispersion parameters. By default is 0.3.
#' @param type type of libaray size, if is "fixed" row sum will be the same or will follows poisson distribution.
#' @return A matrix of n x p each row iid follows Dir(a)
#'
rdm <- function(n, p, lib_size, pi, over_disp = 0.3, type = "fixed"){

  if (missing(lib_size)){
    lib_size = rpois(n, 1e5)
  }

  if (length(lib_size) == 1)
    if(type == "fixed"){
      lib_size <- rep(lib_size, n)
    }else{
      lib_size <- rpois(n, lib_size)
    }

  if (missing(pi))
    pi <- rnorm(p, mean = 14, sd = 4)
  else p <- length(pi)
  pi <- pi/sum(pi)
  P <- rdirichlet(n, pi * (1 - over_disp)/over_disp)
  X <- matrix(0, n, p)

  for (i in 1:n) X[i, ] <- rmultinom(1, lib_size[i], P[i, ])
  list(theta = over_disp, pi = pi, data = X)

}

#' softmax function
#' @param x A vector for compute softmax
#' @return softmax transformation of x
softmax <- function(x){

  return(exp(x)/sum(exp(x)))

}

#' Categorical distribution
#' @param p probability for categorical distribution
#' @return Random number follows Cat(p) distribution
rCate <- function(p){

  return(rmultinom(1, size = 1, prob = p))

}

#' Clipping function to make the value within a range
#' @param x A vector to be clipping
#' @param max_x maximum value of the vector
#' @param min_x minimum value of the vector
#' @return A clipped vector within the range
clipping <- function(x, max_x, min_x){

  x[x > max_x] = max_x

  x[x < min_x] = min_x

  return(x)
}

#' Logit function
#' @param x Input of logit function
#' @return Logistic transformation of the x
logit <- function(x){

  return(1 / (1 + exp(-x)))

}

#' Derivative of logit function
#' @param x A matrix of input for logistic regression
#' @param beta A vector of coefficients.
#' @return \frac{d(logit(x\beta))}{d\beta}
d_logit <- function(x, beta){

  return(x - x * logit(x %*% beta))

}

#' A function for converting a soft probability matrix to factor or one hot matrix.
#' @param prob a matrix to convert.
#' @param idx the type of output. If idx= TRUE will return the index of largest probability of each row. If idx= FALSE will return a one-hot matrix.
#' @return a index of largest probability or one hot matrix.
#' @export

hard_prob <- function(prob, idx = FALSE){

  prob<- as.matrix(prob)
  n <- nrow(prob)
  K <- ncol(prob)

  prob_arg <- apply(prob, 1, which.max)

  if(idx){
    return(prob_arg)
  }else{

    prob_thresh <- matrix(0, nrow = n, ncol = K)

    for(i in 1:n){
      prob_thresh[i,prob_arg[i]] <- 1
    }

    return(prob_thresh)
  }
}

# random sample n samples from a matrix with each row corresponding to the probability of each sample.
#' @param prob a matrix to convert.
#' @return a one-hot matrix follow Cat(prob).
#' @export
#'
sample_prob <- function(prob){

  K <- ncol(prob)

#  Cate_K <- purrr::partial(rCate, K = K)
  sample_res <- t(apply(prob, 1, rCate))
  return(sample_res)
}

#' Inverse function. Using for SAEM optimization.
#' @param it number of iterations
#' @return \frac{1}{it}

inv_func <- function(it){
  return(1/it)
}

#' One hot encoding to convert a factor result to one hot matrix
#' @param x A factor.
#' @return One hot encoding of x.
onehot <- function(x){

  K <- max(x)
  n <- length(x)

  one_hot <- matrix(0, nrow = n, ncol = K)

  for(i in 1:n){
    one_hot[i,x[i]] <- 1
  }

  return(one_hot)
}


#' Transformation of V matrix. From Logistic regression to multinomial regression.
#' @param V Regression coefficient of gating network
#' @param beta_max Maximal value of beta, equivalent to L_inf penalty.
#' @return A coefficient matrix used for multinomial regression.
#'
V_trans <- function(V, beta_max = 10){

  V <- as.matrix(V)
  p <- ncol(V)

  if(p == 1){
    V <- 2 * V
  }else{
    V <- V - V[,1]
  }

  V <- clipping(V, beta_max, -beta_max)

  return(V)
}


#' Transformation of W list. From Logistic regression to Multinomial regression.
#' @param W Regression coefficient of gating network
#' @return A coefficient matrix from W list.
#'
W_trans <- function(W){

  W1 <- list()

  if("matrix" %in% class(W)){

    for(i in 1:ncol(W)){
      W1[[i]] <- W[,i]
    }

  }

  return(W1)
}


#' Transformation of V matrix. From Multinomial regression to Logistic regression.
#' @param beta A vector from logistic regression
#' @param beta_max Maximal value of beta, equivalent to L_inf penalty.
#' @return A coefficient matrix used for logistic regression.

binom_trans <- function(beta, beta_max = 10){

  p <- ncol(beta)

  if(p == 1){
    return(beta)
  }else{

    if(!(all(beta[,1] == 0))){
      beta <- V_trans(beta, beta_max)
    }

    if(p > 2){
      return(beta)
    }else{
      return(beta[,2])
    }

  }
}

#' Function to check if r_i is suitable for input of l layer
#' @param r_i An array for checking if contains any NA values.
#' @param l Number of layers
#' @return An array of (d1, d2, l)
check_r_i_l <- function(r_i, l){

  if(is.na(dim(r_i)[3])){

    r_i0 <- r_i

    d1 <- dim(r_i)[1]
    d2 <- dim(r_i)[2]

    r_i <- array(0, dim = c(d1, d2, l))
    for(i in 1:l){
      r_i[,,i] <- r_i0
    }
  }

  return(r_i)

}

#' Check whether the class of a vector or one hot matrix satisfies a valid class
#' i.e. the number elements within each class larger than min_obs
#' @param y Vector or matrix to check whether the observations of each classes is larger than min_obs
#' @param min_obs Minimal observations for each classes. By default is 7
#' @return A logic value indicate whether all classes in y larger than min_obs.
class_check <- function(y, min_obs = 7){

  y <- as.matrix(y)

  n <- nrow(y)
  K <- ncol(y)

  y <- as.numeric(as.character(y))

  y <- matrix(y, nrow = n, ncol = K)

  if(!length(y)){
    return(FALSE)
  }

  if(K == 1){
    y <- cbind(y, 1 - y)
  }

  y_n <- colSums(y)

  if(any(y_n < min_obs)){
    return(FALSE)
  }else{
    return(TRUE)
  }

}

# Check whether V converge.
.stop_cov <- function(beta_old, beta_new, eps = 1e-6){

  temp1 <- as.matrix(beta_old)
  temp2 <- as.matrix(beta_new)

  return((sum((temp1 - temp2)^2) < eps))

}

# Check whether loss function converge.
.stop_loss <- function(obj, it, eps = 1e-6){

  if(it > 1){

    diff <- obj[it] - obj[it - 1]
    return((diff < eps))

  }else{
    return(FALSE)
  }
}


#' Check whether maximium iteration number reach.
.stop_it <- function(it, itmax){


  return((it >= itmax))

}

# Backtracking of function f.
.backtracking_f <- function(f, x_old, x_new, beta = 0.8, mm = "max"){

  d_x <- x_new - x_old

  f_old <- f(x_old)

  if((is.infinite(f_old)) | (is.nan(f_old))){
    return(x_old)
  }

  it <- 0

  t = 1

  repeat{

    x_new <- x_old + t*d_x

    f_new <- f(x_new)

    if(is.nan(f_new)){
      return(x_old)
    }

    if(mm == "max"){
      cond1 <- ((f_new - f_old) > -1e-9)
    }else{
      cond1 <- ((f_new - f_old) < 1e-9)
    }

    cond2 <- (it > 100)


    if((cond1|cond2)){
      break()
    }

    t <- beta*t
    it <- it + 1
  }

  return(x_new)

}

# Backtracking of W and V simultenousely.

.backtracking_VW <- function(f, W_old, W_new, V_old, V_new, beta = 0.8, mm = "max"){

  K <- length(W_old)

  W_old_df <- as.data.frame(W_old)
  W_new_df <- as.data.frame(W_new)

  p <- nrow(W_new_df)

  d_W_df <- W_new_df - W_old_df
  d_V <- V_new - V_old

  f_old <- f(beta_V = V_old, beta_W = W_old)

  it = 0

  t = 1
  repeat{
    W_new <- as.list(W_old_df + t * d_W_df)
    V_new <- V_old + t*d_V

    f_new <- f(beta_V = V_new, beta_W = W_new)

    if(mm == "max"){
      cond1 <- ((f_new - f_old) > 0)
    }else{
      cond1 <- ((f_new - f_old) < 0)
    }

    cond2 <- (it > 100)

    if(cond1){
      break()
    }
    if(cond2){
      W_new <- W_old
      V_new <- V_old
      break(

      )
    }

    t <- beta*t
    it <- it + 1
  }
  for(i in 1:K){

    W_new[[i]] <- as.matrix(W_new[[i]])

    rownames(W_new[[i]]) <- c("intercept", paste("V",1:(p - 1), sep = ""))

  }
  names(W_new) <- 1:K

  return(list(V_new = V_new, W_new = W_new))
}

# Backtracking of function f.
.backtracking_f_proj <- function(f, x_old, x_new, beta = 0.8, mm = "max",
                                 thresh_cor = thresh_cor, thresh_eig = thresh_eig){

  d_x <- x_new - x_old

  f_old <- f(x_old)

  if((is.infinite(f_old)) | (is.nan(f_old))){
    return(x_old)
  }

  it <- 0

  t = 1

  repeat{

    x_new <- x_old + t*d_x
    x_new <- proj_sd_thresh(x_new, thresh_cor = thresh_cor, thresh_eig = thresh_eig)

    f_new <- f(x_new)

    if(is.nan(f_new)){
      return(x_old)
    }

    if(mm == "max"){
      cond1 <- ((f_new - f_old) > -1e-9)
    }else{
      cond1 <- ((f_new - f_old) < 1e-9)
    }

    cond2 <- (it > 100)


    if((cond1|cond2)){
      break()
    }

    t <- beta*t
    it <- it + 1
  }

  return(x_new)


}

#' hard threshold function
#' @param X Input of hard threshold function
#' @param thresh Threshold
#' @return Hard threshold of X

hard_thresh <- function(X, thresh = 1e-8){

  X[abs(X) < thresh] = 0
  return(X)

}

#' Calculate r_i
.r_i_calc <- function(X, Z, Y, pi_old, W_old, V_old, sgm_old, ratio_max = 20){

  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  nr <- ncol(Y)

  K <- ncol(V_old)

  LL_old <- matrix(0, nrow = n, ncol = K)

  for(i in 1:K){

    LL_old[,i] <- .LLikeli_mvr(X, W_old[,,i], Y = Y, sgm_old[i,], intercept = T, reduce = F)

  }

  pi_ratio <- log(pi_old) - log(pi_old[,1])
  Likeli_ratio <- LL_old - LL_old[,1]

  total_ratio <- pi_ratio + Likeli_ratio

  total_ratio[total_ratio > ratio_max] = ratio_max
  total_ratio[total_ratio < -ratio_max] = -ratio_max

  r_i_new <- t(apply(total_ratio, 1, softmax))

  return(r_i_new)

}

#' Calculate r_i in array form
.r_i_arr_m <- function(X_list, Z, Y_list, pi_old, W_old, V_old, sgm_old, ratio_max){

  l <- length(X_list)

  n <- nrow(Z)
  K <- ncol(V_old)

  r_i_new <- array(1/K,dim = c(n,K,l))

  for(i in 1:l){

    X_temp <- X_list[[i]]
    Y_temp <- Y_list[[i]]
    W_temp <- W_old[[i]]
    sgm_temp <- sgm_old[[i]]

    r_i_new[,,i] <- .r_i_calc(X_temp, Z, Y_temp, pi_old, W_temp, V_old, sgm_temp, ratio_max)

  }
  return(r_i_new)

}

#' Calculate degree of freedom
#' @description This function calculate degree of freedom of NEGMoE object.
#' @param NEMoE_obj A NEGMoE object with fitted result.
#' @param output whether degree of freedom with combined
#'  components or compute within each component.
#' @return a list of degree of freedom in gating network and experts network.
#' @export

calcdf <- function(NEGMoE_obj, output = "all"){

  K <- NEGMoE_obj@K
  gamma <- NEGMoE_obj@NEGMoE_output$gamma
  experts <- NEGMoE_obj@NEGMoE_output$experts

  df_gamma <- .calcdf(gamma, type = "mElnet")
  type <- NEGMoE_obj@type

  if(type == "icov"){

    experts_res <- NEGMoE_obj@NEGMoE_output$experts$Theta

  }else if(type == "mElnet"){

    experts_res <- NEGMoE_obj@NEGMoE_output$experts$W

  }else if(type == "pseudoLL"){

    experts_res <- NEGMoE_obj@NEGMoE_output$experts$W

  }

  df_beta <- matrix(0, nrow = 1, ncol = K)
  for(i in 1:K){

      df_beta[,i] <- .calcdf(experts_res[,,i], type = type)

  }
  if(output == "all"){
    df_beta <- sum(df_beta)
  }

  return(list(df_gamma = df_gamma, df_beta = df_beta))


}

#' Calculate how many variables included in the matrix
#' If version = "all", return how many non-zero rows
#' Else return how many non-zero elements

.calcdf <- function(X, type = "icov", version = "all", eps = 1e-8,
                    intercept = TRUE){


  if(type == "icov"){

    X <- as.matrix(X)
    X_tf <- sum(abs(X[upper.tri(X)]) > eps)
    return(X_tf)

  }else{

    X <- as.matrix(X)

    if(intercept){
      p <- (nrow(X) - 1)
      X <- X[2:(p + 1),,drop = F]
    }else{
      p <- nrow(X)
    }

    if(version == "all"){

      var_tf <- (abs(X) > eps)

      var_tf <- apply(var_tf, 1, any)

      return(sum(as.numeric(var_tf)))

    }else{

      var_tf <- (abs(X) > eps)
      return(sum(apply(var_tf, 2, as.numeric)))
    }
  }



}
