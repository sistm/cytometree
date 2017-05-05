#' TODO
#'@param 
#'
#'@author Chariff Alkhassim
#'
#'@export 
#
Partition2gr <- function(N)
{
  Partitions <- list()
  labelgr1 <- 1
  labelgr2 <- 2
  partition <- rep(labelgr2, N)
  for(i in 1:(N-1))
  {
    temp <- partition
    temp[1:i] <- labelgr1
    Partitions[[i]] <- temp
  }
  Partitions
}
