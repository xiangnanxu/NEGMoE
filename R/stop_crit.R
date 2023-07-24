




# Check whether W converge.

stop_cov_mat <- function(W_old, W_new, eps = 1e-6){

  l <- length(W_old)

  stop_mat <- c()

  for(i in 1:l){
    stop_mat[i] <- stop_cov(W_old[[i]], W_new[[i]])
  }
  return(all(stop_mat))
}

# Check whether loss function converge.

stop_loss <- function(obj, it, eps = 1e-6){

  if(it > 1){

    diff <- obj[it] - obj[it - 1]
    return((diff < eps))

  }else{
    return(FALSE)
  }
}

# Check whether loss function converge, matrix version.

stop_loss_mat <- function(obj, it, eps = 1e-6){

  obj_reduce <- colSums(obj)

  return(stop_loss(obj_reduce, it, eps))
}

# Predict convergence criterion.

stop_pred <- function(obj, it, eps = 0.005){

  if(it > 2){

    c_s <- (obj[it] - obj[it - 1])/(obj[it - 1] - obj[it - 2])
    l_s_tid <- obj[it - 2] + 1/(1 - c_s) * (obj[it - 1] - obj[it - 2])

    return((l_s_tid - obj[it - 1]) < eps)

  }else{
    return(FALSE)
  }

}

