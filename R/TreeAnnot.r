#' Returns the underlying annotation given a tree pattern.
#'
#' @keywords internal
TreeAnnot <- function(labels, combinaisons)
{ 
  if(!is.null(combinaisons))
  {
    out <- c()
    unilabels <- unique(labels)
    for(label in unilabels)
    {
      ind <- labels==label
      ncell <- sum(ind) 
      comb_lab <- combinaisons[match(TRUE, ind),]
      out <- rbind(out, c(comb_lab, label, ncell))  
    }
    out <- out[sort(out[,c(ncol(out))], 
                    decreasing = TRUE, 
                    index.return = TRUE)$ix,]
    out <- cbind(out, round(out[,ncol(out)]/sum(out[,ncol(out)]), 4))
    colnames(out) <- c(colnames(combinaisons), "labels","count","prop")
    as.data.frame(out)
  }
  else
  {
    NULL
  }
}