#' TODO
#'@param 
#'
#'@author Chariff Alkhassim
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
    nn<- length(res$pl_list[[1]])
    for (n in 1:nn)
    {
      minxinf <- min(res$pl_list[[1]][[n]]$x)
      minxsup <- min(res$pl_list[[3]][[n]]$x)  
      maxxinf <- max(res$pl_list[[1]][[n]]$x)
      maxxsup <- max(res$pl_list[[3]][[n]]$y)
      minyinf <- min(res$pl_list[[1]][[n]]$y)
      minysup <- min(res$pl_list[[3]][[n]]$y)
      maxyinf <- max(res$pl_list[[1]][[n]]$y)
      maxysup <- max(res$pl_list[[3]][[n]]$y)
      plot(res$pl_list[[3]][[n]], 
           xlim=c(min(c(minxinf, minxsup)), max(c(maxxinf, maxxsup))),
           ylim=c(min(c(minyinf, minysup)), max(c(maxyinf, maxysup))),
           yaxt = 'n',
           main =res$pl_list[[2]][[n]],
           xlab = res$pl_list[[4]][[n]],
           ylab = "Density",
           col = "blue", lwd=2, type="l")
      lines(res$pl_list[[1]][[n]], col='red', lwd = 2)
    }
  }
  else
  {
    treenodes <- unlist(pl_list[[2]])
    if(sum(nodes%in%treenodes)!=length(nodes))
    {
      wstr <- paste("nodes",paste(c(nodes[as.logical(1-nodes%in%treenodes)]), 
                                  collapse=", "),
                    "are not in the tree", sep = " ")
      stop(wstr)
    }
    inds <- which(treenodes%in%nodes)
    for(ind in inds)
    {
      df <- data.frame(x1 = res$pl_list[[3]][[ind]]$x,
                       y1 = res$pl_list[[3]][[ind]]$y,
                       x2 = res$pl_list[[1]][[ind]]$x,
                       y2 = res$pl_list[[1]][[ind]]$y,
                       Marker = paste0(res$pl_list[[2]][[ind]], 
                                       "\n", res$pl_list[[4]][[ind]]))                           
      p <- ggplot(df, aes(x = x1, y =  y1,col = "1")) +
        geom_line(lwd = 1)+
        geom_line(data = df, aes(x = x2, y =  y2, col = "2"),
                  lwd = 1)+  
        scale_colour_manual(name = "",
                            values = c("blue","red"),
                            labels = c("GM", "KDE"))+ 
        facet_wrap(~ Marker)+
        theme_bw()+
        ylab("Density")+
        xlab("Fluorescence")+
        theme(legend.position="bottom",
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
      print(p)
    }
  }
}





