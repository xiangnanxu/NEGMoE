
# Calculate the number of non-zero columns for a matrix

calc_nz_col <- function(X, eps = 1e-5){
  p <- nrow(X)

  X_tf <- abs(X[2:p,]) > eps
  X_df <- colSums(X_tf)

  return(X_df)
}

# Calculate the number of non-zero rows for a matrix

calc_nz_row <- function(X, eps = 1e-5){
  p <- nrow(X)

  X_tf <- abs(X[2:p,]) > eps
  X_df <- sum(apply(X_tf, 1, any))

  return(X_df)
}

# Calculate the degree of freedom for a MoE result

df_calc <- function(NEGMoE_obj, layered = T, eps = 1e-5){

  experts <- EM_res$V
  W <- EM_res$W

  q <- nrow(V) - 1
  p <- sapply(W, nrow)

  df_V <- calc_nz_row(V)
  if(layered){
    df_W <- sapply(W, calc_nz_row)
  }else{
    W_mat <- list2mat(W)
    df_W <- calc_nz_row(W_mat)
  }

  return(list(df_V = df_V, df_W = df_W))
}

#' Calculate BIC of a RMoE result
#' @param EM_res A list of RMoE result.
#' @param eps A threshold for minimal coefficient for calculate degree of freedom.
#' @return A numerical variable of BIC statistics of the fitted RMoE model.
#' @export

MoE_BIC <- function(EM_res, eps = 1e-5){

  n <- dim(EM_res$r_i)[1]

  V <- EM_res$V
  W <- EM_res$W

  layered = !(is.na(dim(EM_res$r_i)[3]))

  EM_res_df <- df_calc(EM_res, layered)

  if("list" %in% class(EM_res$W)){
    K <- length(EM_res$W)
  }else{
    K = 1
  }

  if(K > 1){
    if(layered){

      LL_temp <- EM_res$obj[,1,]
      LL <- sum(LL_temp[,ncol(LL_temp)])

    }else{
      LL_temp <- EM_res$obj[,1]
      LL <- LL_temp[length(LL_temp)]
    }

    BIC_res <- LL - (EM_res_df$df_V + sum(EM_res_df$df_W)) *(log(n)/2)
  }
  else{

    p <- nrow(EM_res$W)

    LL <- EM_res$LL

    df_W <- sum(abs(EM_res$W[2:p] > eps))

    BIC_res <- LL - df_W*(log(n)/2)
  }

  return(BIC_res)
}


MoE_AIC <- function(EM_res, eps = 1e-5){

  n <- dim(EM_res$r_i)[1]

  V <- EM_res$V
  W <- EM_res$W

  layered = !(is.na(dim(EM_res$r_i)[3]))

  EM_res_df <- df_calc(EM_res, layered)

  if("list" %in% class(EM_res$W)){
    K <- length(EM_res$W)
  }else{
    K = 1
  }

  if(K > 1){
    if(layered){

      LL_temp <- EM_res$obj[,1,]
      LL <- sum(LL_temp[,ncol(LL_temp)])

    }else{
      LL_temp <- EM_res$obj[,1]
      LL <- LL_temp[length(LL_temp)]
    }

    AIC_res <- LL - (EM_res_df$df_V + sum(EM_res_df$df_W))
  }
  else{

    p <- nrow(EM_res$W)

    LL <- EM_res$LL

    df_W <- sum(abs(EM_res$W[2:p] > eps))

    AIC_res <- LL - df_W
  }

  return(AIC_res)
}

#' Calculate ICL of a RMoE result
#' @param EM_res A list of RMoE result.
#' @param eps A threshold for minimal coefficient for calculate degree of freedom.
#' @return A numerical variable of ICL statistics of the fitted RMoE model.
#' @export

MoE_ICL <- function(EM_res, eps = 1e-5){

  n <- dim(EM_res$r_i)[1]

  V <- EM_res$V
  W <- EM_res$W

  layered = !is.na(dim(EM_res$r_i)[3])

  EM_res_df <- df_calc(EM_res, layered)

  if("list" %in% class(EM_res$W)){
    K <- length(EM_res$W)
  }else{
    K = 1
  }

  if(K > 1){
    if(layered){

      LL_temp <- EM_res$obj[,3,]
      LL <- sum(LL_temp[,ncol(LL_temp)])

    }else{
      LL_temp <- EM_res$obj[,3]
      LL <- LL_temp[length(LL_temp)]
    }

    ICL_res <- LL - (EM_res_df$df_V + sum(EM_res_df$df_W)) *(log(n)/2)
  }
  else{

    p <- nrow(EM_res$W)

    LL <- EM_res$LL

    df_W <- sum(abs(EM_res$W[2:p] > eps))

    ICL_res <- LL - df_W*(log(n)/2)
  }


  return(ICL_res)
}

#' Calculate modified ICL of a RMoE result
#' @param EM_res A list of RMoE result.
#' @param eps A threshold for minimal coefficient for calculate degree of freedom.
#' @return A numerical variable of mICL statistics of the fitted RMoE model.
#' @export

MoE_mICL <- function(EM_res, eps = 1e-5){

  n <- dim(EM_res$r_i)[1]

  V <- EM_res$V
  W <- EM_res$W

  layered = !is.na(dim(EM_res$r_i)[3])

  EM_res_df <- df_calc(EM_res, layered)

  if("list" %in% class(EM_res$W)){
    K <- length(EM_res$W)
  }else{
    K = 1
  }

  if(K > 1){

    pi_temp <- colSums(EM_res$pi)
    pi_temp[pi_temp < 1] = 1

    if(layered){

      LL_temp <- EM_res$obj[,3,]
      LL <- sum(LL_temp[,ncol(LL_temp)])

      df_W <- rowSums(sapply(W, calc_nz_col))

    }else{
      LL_temp <- EM_res$obj[,3]
      LL <- LL_temp[length(LL_temp)]

      df_W <- calc_nz_col(list2mat(W))
    }



    mICL_res <- LL - EM_res_df$df_V * log(n)/2 - sum(df_W * log(pi_temp)/2)
  }
  else{
    LL <- EM_res$LL

    p <- nrow(EM_res$W)

    df_W <- sum(abs(EM_res$W[2:p] > eps))

    mICL_res <- LL - df_W*(log(n)/2)
  }

  return(mICL_res)

}
