#' TODO
#'@param 
#'
#'@author Chariff Alkhassim
#'
#'@export 
# 

LeafCenters <- function(CytomeTreeObj)
{
  M <- CytomeTreeObj$M 
  labels <- CytomeTreeObj$labels
  nmar <- ncol(M)
  out <- c()
  Labels <- sort(unique(labels))
  for(mar in 1:nmar)
  {
    mu_leafs <- c()
    for(label in Labels) 
    {
      ind <- which(labels == label)
      mu_leafs <- append(mu_leafs, mean(M[ind, mar]))
    }
    out <- cbind(out, mu_leafs)
  }
  out <- cbind(out, Labels)
  colnames(out) <- c(colnames(M),"leaf")
  out
}
