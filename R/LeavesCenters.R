#' Computes means of leaves given cytomeTreeObj.
#'
#' @keywords internal
LeavesCenters <- function(cytomeTreeObj)
{
  out <- c()
  M <- cytomeTreeObj$M 
  labels <- cytomeTreeObj$labels
  nmar <- ncol(M)
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
