#'Finds the partition which minimize the within-leaves sum of squares.
#'
#'@keywords internal
KmeansOPT<- function(groups, leaves, labels, dat, K)
{
  KK <- length(groups)
  N <- length(labels)
  wlss <- rep(NA, KK)
  if(K == 2)
  {
    for(kk in 1:KK)
    {
      ind1 <- labels%in%leaves[which(groups[[kk]] == 1)]
      ind2 <- as.logical(1-ind1)
      n1 <- sum(ind1)
      n2 <- N-n1
      g1 <- dat[ind1]
      g2 <- dat[ind2]
      d1 <- sum((mean(g1)-g1)**2)
      d2 <- sum((mean(g2)-g2)**2)
      wlss[kk] <- sum(d1, d2)
    }
    winind <- which.min(wlss)
    min_wlss <- wlss[winind]
    return(list("ind" = winind, "val" = min_wlss))
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
      g1 <- dat[ind1]
      g2 <- dat[ind2]
      g3 <- dat[ind3]
      d1 <- sum((mean(g1)-g1)**2)
      d2 <- sum((mean(g2)-g2)**2)
      d3 <- sum((mean(g3)-g3)**2)
      wlss[kk] <- sum(d1, d2, d3)
    }
    winind <- which.min(wlss)
    min_wlss <- wlss[winind]
    return(list("ind" = winind, "val" = min_wlss))
  }
}