#' Plot a heatmap representing the median of the expression for each population and marker using CytomeTree
#' 
#' @param TreeObj An object of class CytomeTree
#' 
#' @param AnnotObj An object of class Annotation
#' 
#' @author Anthony Devaux, Boris Hejblum
#' 
#' @import robustbase pheatmap 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' 
#' @export plot_HMpop
#' 
#' @examples 
#' 
#' # Run CytomeTree
#' data(DLBCL)
#' cellevents <- DLBCL[,c("FL1", "FL2", "FL4")]
#' Tree <- CytomeTree(cellevents, minleaf = 1, t=.1)
#' Annot <- Annotation(Tree,plot=FALSE)
#' 
#' # Plot the cell count
#' plot_HMpop(Tree,Annot)


plot_HMpop <- function(TreeObj, AnnotObj) {
  
  if (class(TreeObj)!="CytomeTree") {
    
    stop("TreeObj must be class of CytomeTree")
    
  }
  
  if (class(AnnotObj)!="Annotation") {
    
    stop("AnnotObj must be class of Annotation")
    
  }
  
  nbClust <- max(AnnotObj$labels)
  data <- matrix(ncol = ncol(TreeObj$M))
  
  for (numclust in 1:nbClust) {
    
    data <- rbind(data,robustbase::colMedians(as.matrix(TreeObj$M[which(AnnotObj$labels==numclust),])))
    
  }
  
  data <- t(data[-1,])
  colnames(data) <- c(1:nbClust)
  
  annotation <- data.frame(Log10_Count = log10(AnnotObj$combinations[,"count"]))
  
  pheatmap::pheatmap(data, scale = "row", annotation_col = annotation, 
           color = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 10, name =
                                                     "RdYlGn"))(100), cluster_rows = FALSE)

}