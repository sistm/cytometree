#' TODO
#'@param 
#'
#'@author Chariff Alkhassim
#'
#'@export 
KmeansOPT<- function(groups, leaves, labels, data, K)
{
  KK <- length(groups)
  N <- length(labels)
  varKK <- rep(NA, KK)
  if(K == 2)
  {
    for(kk in 1:KK)
    {
      ind1 <- labels%in%leaves[which(groups[[kk]] == 1)]
      ind2 <- as.logical(1-ind1)
      n1 <- sum(ind1)
      n2 <- N-n1
      g1 <- data[ind1]
      g2 <- data[ind2]
      d1 <- sum((mean(g1)-g1)**2)
      d2 <- sum((mean(g2)-g2)**2)
      varKK[kk] <- sum(d1, d2) 
    }
    return(list("ind" = which.min(varKK), "val" = min(varKK[kk])))
  }
  else if (K == 3)
  {
    for(kk in 1:KK)
    {
      ind1 <- labels%in%leaves[which(groups[[kk]] == 1)]
      ind2 <- labels%in%leaves[which(groups[[kk]] == 2)]
      ind3 <- labels%in%leaves[which(groups[[kk]] == 3)]
      n1 <- length(ind1)
      n2 <- length(ind2)
      n3 <- length(ind3)
      g1 <- data[ind1]
      g2 <- data[ind2]
      g3 <- data[ind3]
      d1 <- sum((mean(g1)-g1)**2)
      d2 <- sum((mean(g2)-g2)**2)
      d3 <- sum((mean(g3)-g3)**2)
      varKK[kk] <- sum(d1, d2, d3)
    }
    return(list("ind" = which.min(varKK), "val" = min(varKK[kk])))
  }
}