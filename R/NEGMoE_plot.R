###################################################################
# Plot fitting result of gating network
###################################################################
#' Plot gating network
#' @description This function plot the PCA of fitted latent class
#' and their corresponding loadings.
#' @param NEGMoE_obj a NEGMoE object with fitted output.
#' @param PCs Visualization of selected Principal components.
#' @return A graph of PCA plot of nutrition intake and estimated latent classes
#' and its corresponding loadings.
#' @import ggplot2
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_classic coord_flip
#' @export

plotGating <- function(NEGMoE_obj, PCs = c(1,2)){
  if(!length(NEGMoE_obj@NEGMoE_output)){
    error <- "Please fitting NEGMoE before plot."
    return(TRUE)
  }

  gamma <- NEGMoE_obj@NEGMoE_output$gamma

  colnames(gamma) <- paste("latent",1:ncol(gamma), sep = "")
  gamma_df <- as.data.frame(gamma)
  gamma_df <- gamma_df[-1,]
  gamma_df$var_name = rownames(gamma_df)
  gamma_df_melt = reshape2::melt(gamma_df[,-1,drop = F], id.vars = "var_name")
  colnames(gamma_df_melt)[2] = "latent"

  Z <- NEGMoE_obj@Z
  n <- nrow(Z)
  Z1 <- cbind(rep(1,n), Z)
  latent <- t(apply(Z1 %*% gamma,1,softmax))
  latent <- as.factor(hard_prob(latent, idx = TRUE))

  Z_pca <- stats::prcomp(Z, scale. = T)
  Z_pca <- as.data.frame(Z_pca$x[,PCs])
  Z_pca$latent <- latent
  p1 <- ggplot(Z_pca) +
    geom_point(aes(x = Z_pca[,1], y = Z_pca[,2],
                   color = .data$latent))+
    theme_classic() + labs(x = paste0("PC",PCs[1]), y = paste0("PC",PCs[2]))
  p2 <- ggplot(gamma_df_melt) +
    geom_bar(aes(x = .data$var_name, y = .data$value, fill = .data$latent),
             stat = "identity") + coord_flip() +
    labs(x= "", y ="") + theme_classic()

  return(list(p1, p2))

}
###################################################################
# Plot fitting result of experts network
###################################################################
#' Plot experts network
#' @description This function plot the estimated coefficients of each latent class
#' on each level.
#' @param NEGMoE_obj a NEGMoE object with fitted output.
#' @return A list of graph of fitted coefficients of each latent class
#' on each level.
#' @import ggplot2
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_classic
#' @import ggnetwork

#' @export

plotExperts <- function(NEGMoE_obj, thresh = 1e-2, pseudo_type = "or"){


  type <- NEGMoE_obj@type
  K <- NEGMoE_obj@K

  experts_res <- NEGMoE_obj@NEGMoE_output$experts
  nodes <- data.frame(name = "")
  edge_list <- list()

  if(type == "icov"){
    Theta <- experts_res$Theta
    for(k in 1:K){
      Theta_01 <- matrix(0, nrow = nrow(Theta[,,k]), ncol = ncol(Theta[,,k]))
      Theta_01[abs(Theta[,,k]) > thresh] = sign(Theta[,,k][abs(Theta[,,k]) > thresh])
      Theta_01 <- (Theta_01 + t(Theta_01))/2
      rownames(Theta_01) <- rownames(Theta[,,k])
      colnames(Theta_01) <- colnames(Theta[,,k])
      Theta_01[lower.tri(Theta_01,diag = TRUE)] = 0
      Theta_tmp <- reshape2::melt(Theta_01)
      colnames(Theta_tmp) <- c("from", "to", "type")
      Theta_tmp$type <- as.factor(Theta_tmp$type)
      Theta_tmp$from <- as.character(make.names(Theta_tmp$from))
      Theta_tmp$to <- as.character(make.names(Theta_tmp$to))
      Theta_tmp$latent = k
      edge_list[[k]] <- Theta_tmp
      nodes_list <- c(Theta_tmp$from, Theta_tmp$to)
      nodes <- rbind(nodes, data.frame(name = unique(nodes_list)))
    }
  }else{
    W <- experts_res$W[-1,,]

    for(k in 1:K){
      W_01 <- matrix(0, nrow = nrow(W[,,k]), ncol = ncol(W[,,k]))
      W_01[abs(W[,,k]) > thresh] = sign(W[,,k][abs(W[,,k]) > thresh])
      W_01 <- (W_01 + t(W_01))/2
      rownames(W_01) <- rownames(W[,,k])
      colnames(W_01) <- colnames(W[,,k])
      W_01[lower.tri(W_01,diag = TRUE)] = 0
      W_tmp <- reshape2::melt(W_01)
      if(pseudo_type == "or"){
        W_tmp = W_tmp[W_tmp$value != 0,]
        W_tmp$value = sign(W_tmp$value)
      }else{
        W_tmp = W_tmp[abs(W_tmp$value) == 1,]
        W_tmp$value = sign(W_tmp$value)
      }
      colnames(W_tmp) <- c("from", "to", "type")
      W_tmp$type <- as.factor(W_tmp$type)
      W_tmp$from <- as.character(make.names(W_tmp$from))
      W_tmp$to <- as.character(make.names(W_tmp$to))
      W_tmp$latent = k
      edge_list[[k]] <- W_tmp
      nodes_list <- c(W_tmp$from, W_tmp$to)
      nodes <- rbind(nodes, data.frame(name = unique(nodes_list)))
    }
  }
  nodes <- nodes[-1,,drop = F]
  nodes <- nodes[!duplicated(nodes$name),,drop = F]
  nodes$name <- as.character(nodes$name)
  nodes$label <- nodes$name

#  gen_id <- match(nodes$name, make.names(taxa_table@.Data[,6]))
#  nodes$phy <- unname(taxa_table@.Data[gen_id, 2])

  edge_list <- do.call(rbind,edge_list)
  net <- network::network(edge_list, matrix.type = "edgelist", directed = FALSE,
                          vertices = nodes, multiple = TRUE)

  p <- ggplot(
    ggnetwork(net, by = "latent"),
    aes(x = x, y = y, xend = xend, yend = yend)
  ) +
    geom_edges(aes(linetype = type), size = 1,
               alpha=0.5) +
    geom_nodes(aes(color = phy), size = 5, alpha = 0.5)  +
    geom_nodetext(aes(label = label), size = 3) +
    scale_linetype_manual(values = c("1" = "solid", "-1" = "dashed")) +
    #  geom_edgetext(aes(label = round(log.p.value,2)), color = "grey25") +
    theme_blank() + guides() + facet_wrap(~latent)
  return(p)
}

###################################################################
# Plot fitting result of top paired differential correlation
###################################################################
#' Plot experts network
#' @description This function do scatter plot for differential correlations.
#' @param NEGMoE_obj a NEGMoE object with fitted output.
#' @return A list of graph of fitted coefficients of each latent class
#' on each level.
#' @import ggplot2
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_classic
#' @export
#'
plotScatter <- function(NEGMoE_obj, topdiff = 3, Varlist = NULL){

  type <- NEGMoE_obj@type
  K <- NEGMoE_obj@K

  experts_res <- NEGMoE_obj@NEGMoE_output$experts
  if(type == "icov"){

    Theta <- experts_res$Theta
    Normalized_Theta <- Theta
    for(k in 1:K){
      Normalized_Theta[,,k] <- stats::cov2cor(Theta[,,k])
    }
    Normalized_var <- apply(Normalized_Theta, c(1,2), var)


  }else{

    W <- experts_res$W
    Normalized_W <- W[-1,,]
    for(k in 1:K){
      Normalized_W[,,k] <- Normalized_W[,,k]
    }
    Normalized_var <- apply(Normalized_W, c(1,2), var)
  }
  Normalized_var <- (Normalized_var + t(Normalized_var))/2
  Normalized_var[lower.tri(Normalized_var)] <- 0
  Normalized_var_melt <- reshape2::melt(Normalized_var)
  Normalized_var_melt <- Normalized_var_melt[Normalized_var_melt$value != 0,]
  Normalized_var_melt <- Normalized_var_melt[order(Normalized_var_melt$value,
                                                   decreasing = T),]
  latent <- as.factor(hard_prob(NEGMoE_obj@NEGMoE_output$pi, idx = TRUE))
  p_list <- list()
  cor_df <- data.frame(Var1 = "", Var2 = "", value = 0, latent = 0, method = "")
  if(!is.null(Varlist)){
    topdiff = nrow(Varlist)
    Normalized_var_melt = Varlist
  }
  for(i in 1:topdiff){
    df <- data.frame(x = NEGMoE_obj@X[,Normalized_var_melt[i,1]],
                     y = NEGMoE_obj@X[,Normalized_var_melt[i,2]],
                     latent = latent)
    p_list[[i]] <- ggplot(df, aes(x = .data$x, y = .data$y,
                        color = .data$latent)) + geom_point() +
                  facet_wrap(~.data$latent, scales = "free") +
                  geom_smooth(method = "lm") +
                  labs(x = Normalized_var_melt[i,1], y = Normalized_var_melt[i,2]) +
                  theme_bw()
    cor_tmp <- data.frame(Var1 = Normalized_var_melt[i,1], Var2 = Normalized_var_melt[i,2],
                          value = c(cor(df$x, df$y, method = "pearson"),
                                    cor(df$x, df$y, method = "spearman")),
                          latent = "all",
                          method = c("pearson","spearman"))
    cor_df <- rbind(cor_df, cor_tmp)
    for(k in 1:K){
      cor_tmp <- data.frame(Var1 = Normalized_var_melt[i,1], Var2 = Normalized_var_melt[i,2],
                            value = c(cor(df$x[df$latent == k], df$y[df$latent == k], method = "pearson"),
                                      cor(df$x[df$latent == k], df$y[df$latent == k], method = "spearman")),
                            latent = k,
                            method = c("pearson","spearman"))
      cor_df <- rbind(cor_df, cor_tmp)
    }
  }
  cor_df <- cor_df[-1,]


  return(list(plots = p_list, cor = cor_df))
}
