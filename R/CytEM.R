CytEM <- function(M, indices, minleaf, level, t)
{
  if(class(M)!="matrix")
  {
    M <- as.matrix(M)
  }
  n <- nrow(M)
  p <- ncol(M)
  if(n <= minleaf)
  {
    return(list("mark_not_dis" = 1:p))
  }
  mu1 <- mu2 <- Var1 <- Var2 <- pi1 <- pi2 <- c()
  Resaic <- ind_marker <- mark_not_dis <- c()
  child <- list()
  for(j in 1:p)
  {
    flag_uni <- 0
    M_j <- M[,j]
    mc_uni <- Mclust(M_j, 1)
    mc_mix <- Mclust(M_j, 2, modelNames = "E")
    ind1 <- which(mc_mix$classification == 1)
    ind2 <- which(mc_mix$classification == 2)
    m_zeros <- 0
    thresh_zeros <- .25   # empirically determined
    if(length(which(!M_j)) > length(M_j) * thresh_zeros)
    {
      m_zeros <- 1
    }
    if(length(ind1) < minleaf | length(ind2) < minleaf | m_zeros)
    {
      mark_not_dis <- append(mark_not_dis, j)
      next()
    }
    M1 <- M_j[ind1]
    M2 <- M_j[ind2]
    aic_uni <- 2*mc_uni$df - 2*mc_uni$loglik
    aic_mix <- 2*mc_mix$df - 2*mc_mix$loglik
    aic_norm <- (aic_uni - aic_mix)/n
    flagComp <- 0      
    if(flag_uni | aic_norm < t)
    {
      mark_not_dis <- append(mark_not_dis, j)
    }
    else
    {
      label <- mc_mix$classification
      mean_M1 <- mean(M1)
      mean_M2 <- mean(M2)
      var_M1 <- var(M1)
      var_M2 <- var(M2)
      pi_M1 <- length(M1)/n
      pi_M2 <- 1 - pi_M1 
      ind_marker <- append(ind_marker,j)
      if(mean_M1 > mean_M2)
      {
        label[ind1] <- 1
        label[ind2] <- 0
        mu1 <- append(mu1, mean_M2)
        mu2 <- append(mu2, mean_M1)
        Var1 <- append(Var1, var_M2)
        Var2 <- append(Var2, var_M1)
        pi1 <- append(pi1, pi_M2)
        pi2 <- append(pi2, pi_M1)
      }
      else
      {
        label[ind1] <- 0
        label[ind2] <- 1
        mu1 <- append(mu1, mean_M1)
        mu2 <- append(mu2, mean_M2)
        Var1 <- append(Var1, var_M1)
        Var2 <- append(Var2, var_M2)
        pi1 <- append(pi1, pi_M1)
        pi2 <- append(pi2, pi_M2)
        
      }
      lve <- length(Resaic)
      if(!lve)
      {
        child$L <-  indices[which(label == 0)]
        child$R <-  indices[which(label == 1)]
      } 
      else 
      {
        if(max(Resaic) < aic_norm)
        {
          child$L <-  indices[which(label == 0)]
          child$R <-  indices[which(label == 1)]
        }
      }
      Resaic <-  append(Resaic, aic_norm)
    } 
  }
  len_ind_marker <- length(ind_marker)
  if (!len_ind_marker)
  {
    return(list("mark_not_dis"=mark_not_dis))
  }
  else if(len_ind_marker == 1)
  {
    return(list("nAIC" = Resaic, "ind" = ind_marker, 
                "mark_not_dis" = mark_not_dis, "child" = child, 
                "mu1" = mu1, "mu2" = mu2, 
                "Var1" = Var1, "Var2" = Var2, 
                "pi1" = pi1, "pi2" = pi2))
  }
  else
  {
    res <- cbind(Resaic, ind_marker, mu1, mu2, Var1, Var2, pi1, pi2)
    res <- res[order(res[,1],decreasing = TRUE),]
    return(list("nAIC" = res[,1], "ind" = res[,2], 
                "mark_not_dis" = mark_not_dis, "child" = child, 
                "mu1"= res[1,3], "mu2"= res[1,4],
                "Var1" = res[1,5], "Var2" = res[1,6],
                "pi1"= res[1,7], "pi2" = res[1,8]))
  }
}



