clusTag <- function(labels, combinations, Annotate)
{  
  out <- c()
  unilabels <- sort(unique(labels))
  for(label in unilabels)
  {
    ind <- labels==label
    ncell <- sum(ind) 
    comb_lab <- combinations[match(TRUE, ind),]
    out <- rbind(out, c(comb_lab, label, ncell))  
  }
  out <- out[sort(out[,c(ncol(out))], 
                  decreasing = TRUE, 
                  index.return = TRUE)$ix,]
  out <- cbind(out, round(out[,ncol(out)]/sum(out[,ncol(out)]), 4))
  colnames(out) <- c(colnames(combinations), "labels","ncells","prop")
  out
}



