#'E-M algorithm
#'
#'@import mclust
#'@importFrom stats var
#'@keywords internal

CytEM <- function(M, indices, minleaf, level, t)
{
  if(class(M)!="matrix")
  {
    M <- as.matrix(M)
  }
  n <- nrow(M)
  p <- ncol(M)
  if(n <= minleaf )
  {
    return(list("mark_not_dis" = 1:p))
  }
  nEMdegenerate <- 5
  if(minleaf < nEMdegenerate)
  {
    minleaf <- nEMdegenerate
  }
  parameters <- aic_norm_old <- mark_not_dis <- c()
  t_zer <- .25
  child <- list()
  for(j in 1:p)
  {
    flag_uni <- 0
    M_j <- M[,j]
    mc_uni <- Mclust(M_j, 1)
    mc_mix <- Mclust(M_j, 2, modelNames = "E")
    ind1 <- which(mc_mix$classification == 1)
    ind2 <- which(mc_mix$classification == 2)
    if(length(ind1)<minleaf|length(ind2)<minleaf|is.null(mc_mix)
       |length(which(!M_j))>(n*t_zer))
    {
      mark_not_dis <- append(mark_not_dis, j)
      next()
    }
    M1 <- M_j[ind1]
    M2 <- M_j[ind2]
    aic_uni <- 2*mc_uni$df - 2*mc_uni$loglik
    aic_mix <- 2*mc_mix$df - 2*mc_mix$loglik
    aic_norm_new <- (aic_uni - aic_mix)/n
    flagComp <- 0
    if(flag_uni | aic_norm_new < t)
    {
      mark_not_dis <- append(mark_not_dis, j)
    }
    else
    {
      label <- mc_mix$classification
      mean_M1 <- mean(M1)
      mean_M2 <- mean(M2)
      pi_M1 <- length(M1)/n
      pi_M2 <- 1 - pi_M1
      if(mean_M1 > mean_M2)
      {
        label[ind1] <- 1
        label[ind2] <- 0
        temparameters <- c(aic_norm_new, j, mean_M2, mean_M1,
                           stats::var(M2), stats::var(M1), pi_M2, pi_M1)
      }
      else
      {
        label[ind1] <- 0
        label[ind2] <- 1
        temparameters <- c(aic_norm_new, j, mean_M1, mean_M2,
                           stats::var(M1), stats::var(M2), pi_M1, pi_M2)
      }
      if(is.null(aic_norm_old))
      {
        child$L <-  indices[which(label == 0)]
        child$R <-  indices[which(label == 1)]
        aic_norm_old <- aic_norm_new
      }
      else
      {
        if(aic_norm_old < aic_norm_new)
        {
          child$L <-  indices[which(label == 0)]
          child$R <-  indices[which(label == 1)]
          aic_norm_old <- aic_norm_new
        }
      }
      parameters <- rbind(parameters, temparameters)
    }
  }
  nnrowpara <- nrow(parameters)
  if(is.null(nnrowpara))
  {
    return(list("mark_not_dis"=mark_not_dis))
  }
  if(nnrowpara > 1)
  {
    parameters <- parameters[order(parameters[,1],decreasing = TRUE),]
  }
  return(list("mark_not_dis" = mark_not_dis, "child" = child,
              "nAIC" = parameters[,1], "ind" = parameters[,2],
              "mu1"= parameters[1,3], "mu2"= parameters[1,4],
              "Var1" = parameters[1,5], "Var2" = parameters[1,6],
              "pi1"= parameters[1,7], "pi2" = parameters[1,8]))
}

