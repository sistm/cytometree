#' Plot the cell count for each population using CytomeTree.
#' 
#' @param AnnotObj An object of class Annotation.
#' 
#' @param nbpop Number indicating the maximum of population plotted.
#' Defaut is \code{10}
#' 
#' @param mincount Number indicating the minimum of cell count
#' for the populations. Defaut is \code{1}.
#' 
#' @param maxcount Number indicating the maximum of cell count
#' for the populations. Defaut is \code{NULL} i.e no maximum selected.
#' 
#' @param y_axis a character string either \code{"abs_count"} or \code{"prop"} indicating 
#' whether the absolute cell count or the relative populations proportions should be ploted.
#' Default is \code{"abs_count"}.
#' 
#' @author Anthony Devaux, Boris Hejblum
#' 
#' @import ggplot2
#' 
#' @export
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
#' plot_cytopop(Annot)


plot_cytopop <- function(AnnotObj, nbpop = 10, mincount = 1, maxcount = NULL, y_axis=c("abs_count", "prop")) {
  
  if(length(y_axis)>1){
    y_axis <- y_axis[1]
  }

  if (class(AnnotObj)!="Annotation") {
    
    stop("AnnotObj must be class of Annotation")
    
  }
  
  if (!is.null(nbpop)) {
    
    if (class(nbpop)!="numeric") {
      
      stop("nbpop must be class of numeric")
      
    }
    
  }
  
  if (class(mincount)!="numeric") {
    
    stop("mincount must be class of numeric")
    
  }
  
  if (!is.null(maxcount)) {
    
    if (class(maxcount)!="numeric") {
      
      stop("maxcount must be class of numeric")
      
    }else{
      
      if (maxcount<=mincount) {
        
        stop("maxcount must be higher than mincount")
        
      }
      
    }
    
  }
    
  data <- data.frame(AnnotObj$combinations[,c("leaves","count", "prop")])
  data <- subset(data, data[,"count"] >= mincount)
  
  if (!is.null(maxcount)) {
    
    data <- subset(data, data[,"count"] < maxcount)
    
  }
  
  
  if(y_axis == "abs_count"){
    data <- subset(data, select=c("leaves", "count"))
  }
  else if(y_axis == "prop"){
    data <- subset(data, select=c("leaves", "prop"))
  }
  
  
  if (!is.null(nbpop)) {
    
    if (nbpop<dim(data)[1]) {
      
      data <- data[1:nbpop,]
      
    }
    
  }
  
  if (dim(data)[1]!=0) {
    
    p <- ggplot(data = data)
    
    if(y_axis == "abs_count"){
      p <- p + 
        geom_bar(aes(x=as.factor(leaves), y=count), stat = "identity", fill = "steelblue") +
        scale_y_log10() +
        ylab("Cell count")
    }
    else if(y_axis == "prop"){
      p <- p + 
        geom_bar(aes(x=as.factor(leaves), y=prop), stat = "identity", fill = "steelblue") +
        ylab("Poulation proportion") 
    }
    
    p <- p +
      scale_x_discrete(limits = factor(data$leaves)) +
      xlab("Populations") +
      theme(axis.title=element_text(size=15),
            axis.text=element_text(size=12,face = "bold")) +
      theme_bw()
    
    
    print(p)
    
  
  }else{
    
    stop("No population found... Consider using other settings")
    
  }

}