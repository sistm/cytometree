#' Computes (N-1)*(N-2) possible partitions of N into 3 subsets. 
#'
#' @keywords internal
Partition3gr <- function(N)
{
  partitions <- list()
  partition <- rep(NA, N)
  loop1 <- 1:(N-2)
  cpt <- 1
  labelgr1 <- 1
  labelgr2 <- 2
  labelgr3 <- 3
  for(i in loop1)
  {
    partition[1:i] <- labelgr1
    loop2 <- (i+1):(N-1) 
    for(j in loop2)
    {
      partition[(i+1):j] <- labelgr2
      partition[(j+1):N] <- labelgr3
      partitions[[cpt]] <- partition
      cpt <- cpt + 1
      
    }
  }
  partitions
}




