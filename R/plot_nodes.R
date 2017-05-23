#' Plot the distribution of the observed cells at each node 
#' of the binary tree built using CytomeTree.
#' 
#'@param CytomeTreeObj An object of class CytomeTree.
#'
#'@param nodes A vector of class character containing the name of
#'nodes for which the distribution is to be plotted. Default is 
#'\code{NULL}, and plots the distribution of each node.  
#'
#'@author Chariff Alkhassim
#'
#'@import ggplot2 graphics 
#'
#'@export
plot_nodes <- function(CytomeTreeObj, nodes=NULL)
{
  if(class(CytomeTreeObj) != "CytomeTree")
  {
    stop("CytomeTreeObj must be of class CytomeTree")
  }
  if(!is.null(nodes))
  {
    if(class(nodes)!="character")
    {
      stop("nodes must be of class character")
    }
  }
  pl_list <- CytomeTreeObj$pl_list
  if(is.null(nodes))
  {
    nn<- length(pl_list[[1]])
    for (n in 1:nn)
    {
      minxinf <- min(pl_list[[1]][[n]]$x)
      minxsup <- min(pl_list[[3]][[n]]$x)  
      maxxinf <- max(pl_list[[1]][[n]]$x)
      maxxsup <- max(pl_list[[3]][[n]]$y)
      minyinf <- min(pl_list[[1]][[n]]$y)
      minysup <- min(pl_list[[3]][[n]]$y)
      maxyinf <- max(pl_list[[1]][[n]]$y)
      maxysup <- max(pl_list[[3]][[n]]$y)
      graphics::plot(pl_list[[3]][[n]], 
                     xlim=c(min(c(minxinf, minxsup)), max(c(maxxinf, maxxsup))),
                     ylim=c(min(c(minyinf, minysup)), max(c(maxyinf, maxysup))),
                     yaxt = 'n',
                     main =pl_list[[2]][[n]],
                     xlab = pl_list[[4]][[n]],
                     ylab = "Density",
                     col = "blue", lwd=2, type="l")
      lines(pl_list[[1]][[n]], col='red', lwd = 2)
    }
  }
  else
  {
    treenodes <- unlist(pl_list[[2]])
    if(sum(nodes%in%treenodes)!=length(nodes))
    {
      logicalind <- as.logical(1-nodes%in%treenodes)
      if(length(logicalind) >1)
      {
        wstr <- paste("nodes",paste(c(nodes[]), 
                                    collapse=", "),
                      "are not in the tree", sep = " ")
      }
      else
      {
        wstr <- paste("node",paste(c(nodes[]), 
                                   collapse=", "),
                      "is not in the tree", sep = " ")
      }
      stop(wstr)
    }
    inds <- which(treenodes%in%nodes)
    for(ind in inds)
    {
      df <- data.frame("x1" = pl_list[[3]][[ind]]$x,
                       "y1" = pl_list[[3]][[ind]]$y,
                       "x2" = pl_list[[1]][[ind]]$x,
                       "y2" = pl_list[[1]][[ind]]$y,
                       "Marker" = paste0(pl_list[[2]][[ind]], 
                                         "\n", pl_list[[4]][[ind]]),
                       "colx1" = as.factor(rep(1,length(pl_list[[3]][[ind]]))),
                       "colx2" = as.factor(rep(2,length(pl_list[[3]][[ind]]))))  
      p <- ggplot(df, ggplot2::aes_string(x = "x1", y = "y1",colour = "colx1"))+
        ggplot2::geom_line(lwd = 1)+
        ggplot2::geom_line(data = df, ggplot2::aes_string(x = "x2", y = "y2", 
                                                          colour = "colx2"),
                           lwd = 1)+  
        ggplot2::scale_colour_manual(name = "",
                                     values = c("blue","red"),
                                     labels = c("GM", "KDE"))+ 
        ggplot2::facet_wrap(~ Marker)+
        ggplot2::theme_bw()+
        ggplot2::ylab("Density")+
        ggplot2::xlab("Fluorescence")+
        ggplot2::theme(legend.position="bottom",
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())
      print(p)
    }
  }
}





