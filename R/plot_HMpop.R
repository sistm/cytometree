#' Plot a heatmap representing the median of the expression for each population and marker using CytomeTree
#' 
#' Under construction...
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
#' @importFrom circlize colorRamp2
#' 
#' @export plot_HMpop
#' 
#' @examples 
#' \dontrun{
#' # Run CytomeTree
#' data(DLBCL)
#' cellevents <- DLBCL[,c("FL1", "FL2", "FL4")]
#' Tree <- CytomeTree(cellevents, minleaf = 1, t=.1)
#' Annot <- Annotation(Tree,plot=FALSE)
#' 
#' # Plot the cell count
#' plot_HMpop(Tree,Annot)
#' }


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
  
  color = colorRamp2(seq(min(data), max(data), length = 3), c("red", "#EEEEEE", "green"), space = "RGB")
  
  df <- data.frame(Log10_Count = log10(AnnotObj$combinations[,"count"]))
  
  annotation <- HeatmapAnnotation(df = df,
                                  col = list(Log10_Count = circlize::colorRamp2(c(min(df), max(df)), 
                                                                      c("yellow", "orange")))) 
  
  Heatmap(data, col = color, name = "Median of intensity by population", 
          column_title = "Populations", row_title = "Markers",
          column_title_side = "bottom",
          top_annotation = annotation,
          cluster_rows = FALSE,
          row_names_side = "left",
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.1f", data[i, j]), x, y, gp = gpar(fontsize = 10))
          })
  
  
  
  ################################################################################
  
  #annotation <- data.frame(Log10_Count = log10(AnnotObj$combinations[,"count"]))
  
  #pheatmap::pheatmap(data, scale = "row", annotation_col = annotation, 
  #         color = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 10, name =
  #                                                   "RdYlGn"))(100), cluster_rows = FALSE)

}