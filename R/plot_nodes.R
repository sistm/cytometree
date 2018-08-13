#' Plot the distribution of the observed cells at each node 
#' of the binary tree built using CytomeTree.
#' 
#'@param CytomeTreeObj An object of class CytomeTree.
#'
#'@param nodes A vector of class character containing the name of
#'the nodes for which the distribution is to be plotted. Default is 
#'\code{NULL}, and plots the distribution for each node.  
#'
#'@param nodesPerCol an integer specifying the number of plots to be
#'displayed per column when plotting multiple nodes at once. Default is 
#'\code{NULL}.
#'
#'@param nodesPerRow an integer specifying the number of plots to be
#'displayed per row when plotting multiple nodes at once. Default is 
#'\code{NULL}.
#'
#'@param ... further arguments to be passed to \code{\link[cowplot]{plotgrid}}.
#'
#'@details if both \code{nodesPerCol} and \code{nodesPerRow} are \code{NULL}
#'then all the nodes are plotted on a single page.
#'
#'@return a list of \code{ggplot2} plot objects, containing each node plot.
#'
#'@author Chariff Alkhassim, Boris Hejblum
#'
#'@examples
#'
#'data(DLBCL)
#'myct <- CytomeTree(DLBCL[, c("FL1", "FL2", "FL4")], minleaf = 1, t=.1)
#'plot_nodes(myct)
#'
#'@import ggplot2 graphics 
#'@importFrom cowplot plot_grid
#'
#'@export
plot_nodes <- function(CytomeTreeObj, nodes=NULL, nodesPerCol = NULL, 
                       nodesPerRow = NULL, ...){
  
  if(class(CytomeTreeObj) != "CytomeTree"){
    stop("'CytomeTreeObj' must be of class CytomeTree")
  }
  if(!is.null(nodes)){
    if(class(nodes) != "list"){
      stop("'nodes' argument must be a list")
    }else{
      if(class(unlist(nodes)) != "character")
        stop("elements 'nodes' must be of class character")
    }
  }
  
  pl_list <- CytomeTreeObj$pl_list
  if(is.null(nodes)){
    nodes <- unlist(pl_list$nodes)
  }
  nNodes <- length(nodes)
  
  
  treenodes <- unlist(pl_list$nodes)
  if(sum(nodes%in%treenodes)!=length(nodes)){
    logicalind <- as.logical(1-nodes%in%treenodes)
    if(sum(logicalind) >1){
      wstr <- paste("Nodes", paste(c(nodes[logicalind]), collapse=", "),
                    "are not in 'CytomeTreeObj'", sep = " ")
    }else{
      wstr <- paste("Node", paste(c(nodes[logicalind]), collapse=", "),
                    "is not in 'CytomeTreeObj'", sep = " ")
    }
    stop(wstr)
  }
  
  inds <- which(treenodes%in%nodes)
  df <- list()
  plot_list <- list()
  for(ind in inds){
    df[[ind]] <- data.frame("x" = c(pl_list$gmm[[ind]]$x, 
                                    pl_list$kde[[ind]]$x),
                            "y" = c(pl_list$gmm[[ind]]$y,
                                    pl_list$kde[[ind]]$y),
                            "Marker" = rep(paste0(pl_list$nodes[[ind]], "\n", pl_list$legend[[ind]]),
                                           length(pl_list$gmm[[ind]]$x) + length(pl_list$kde[[ind]]$x)),
                            "Estimator" = factor(c(rep("GMM", length(pl_list$gmm[[ind]]$x)),
                                                   rep("KDE", length(pl_list$kde[[ind]]$x)))
                            )
    )
    p <- ggplot(df[[ind]], ggplot2::aes_string(x = "x", y = "y")) +
      ggplot2::geom_line(ggplot2::aes_string(colour = "Estimator"), lwd = 1) +  
      ggplot2::geom_polygon(ggplot2::aes_string(fill = "Estimator"), alpha=0.3) +
      ggplot2::scale_colour_manual(name = "",
                                   values = c("blue","red"),
                                   labels = c("GM", "KDE")) + 
      ggplot2::scale_fill_manual(name = "",
                                 values = c("lightblue","darkred"),
                                 labels = c("GM", "KDE")) + 
      ggplot2::facet_wrap(~ Marker, nrow = nodesPerRow, ncol = nodesPerCol, 
                          scales = "free", drop=FALSE) +
      ggplot2::theme_bw() +
      ggplot2::ylab("Density") +
      ggplot2::xlab("Fluorescence") +
      ggplot2::theme(legend.position="bottom",
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    
    plot_list[[ind]] <- p
  }
  
  browser()
  if(length(nodes)>1 & (!is.null(nodesPerRow) | !is.null(nodesPerCol))){
    if(is.null(nodesPerRow)){
      nodesPerRow <- 1
    }else if(is.null(nodesPerCol)){
      nodesPerCol <- 1
    }
    nodesPerPage <- nodesPerRow * nodesPerCol
    nbPages <- ceiling(length(nodes)/nodesPerPage)
    
    if(nbPages >1 & nbPages != nNodes/nodesPerPage){
      for(i in 1:(nbPages)){
        print(cowplot::plot_grid(plotlist = plot_list[((i-1)*nodesPerPage + 1):min(nNodes, (i*nodesPerPage))], 
                           nrow = nodesPerCol, ncol = nodesPerRow))
      }
    }else{
      print(cowplot::plot_grid(plotlist = plot_list, nrow = nodesPerRow, ncol = nodesPerCol, ...))
    }
  }else{
    print(cowplot::plot_grid(plotlist = plot_list, nrow = nodesPerRow, ncol = nodesPerCol, ...))
  }
  invisible(plot_list)
}
