#' Plot the cell count for each population using CytomeTree.
#' 
#' @param AnnotObj An object of class Annotation.
#' 
#' @param nbpop Number indicating the number of population plotted.
#' Defaut is \code{10}
#' 
#' @param mincount Number indicating the minimum of cell count
#' for the populations. Defaut is \code{1}.
#' 
#' @param maxcount Number indicating the maximum of cell count
#' for the populations. Defaut is \code{NULL} i.e no maximum selected.
#' 
#' @author Anthony DEVAUX
#' 
#' @import ggplot2
#' 
#' @export


plot_cytopop <- function(AnnotObj, nbpop = 10, mincount = 1, maxcount = NULL) {

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
    
  data <- data.frame(AnnotObj$combinations[,c("leaves","count")])
  data <- subset(data, data[,"count"] >= mincount)
  
  if (!is.null(maxcount)) {
    
    data <- subset(data, data[,"count"] < maxcount)
    
  }
  
  if (!is.null(nbpop)) {
    
    data <- data[1:nbpop,]
    
  }
  
  if (dim(data)[1]!=0) {
  
    p <- ggplot(data = data) +
      geom_bar(aes(x=as.factor(leaves),y=count), stat = "identity", fill = "steelblue") +
      scale_y_log10() +
      scale_x_discrete(limits = factor(data$leaves)) +
      xlab("Populations") +
      ylab("Cell count") +
      theme(axis.title=element_text(size=15),
            axis.text=element_text(size=12,face = "bold"))
    
    print(p)
  
  }else{
    
    stop("No population found ! Change settings please.")
    
  }

}