#' Plot the binary tree built using CytomeTree.
#' 
#'@param CytomeTreeObj An object of class CytomeTree.
#'
#'@param Ecex Number indicating the amount by which text
#'on the edges should be scaled. Default is \code{1}.
#'
#'@param Ecolor An integer or a string of character
#'to color edges of the graph.  Default is \code{8}.
#'
#'@param Vcex Number indicating the amount by which text
#'in the vertices should be scaled.  Default is \code{.8}.  
#'
#'@param Vcolor A vector of class numeric or character to color
#'vertices of the graph. Default is \code{0}.
#'
#'@param ... additional arguments to be passed to \code{\link{plot_graph}}
#'
#'@author Chariff Alkhassim
#'
#'@import igraph
#'
#'@export 
plot_graph <- function(CytomeTreeObj, Ecex = 1, Ecolor = 8,
                       Vcex = .8, Vcolor = 0, ...)
{
  if(class(CytomeTreeObj) != "CytomeTree")
  {
    stop("CytomeTreeObj must be of class CytomeTree.")
  }
  Tree <- CytomeTreeObj$mark_tree
  if(is.null(Tree))
  {
    return(cat("CytomeTree found a single population.\n"))
  }
  Tree_level <- length(Tree)
  adj_list <- c()
  cptleaf <- 1
  for(level in 1:(Tree_level - 1))
  {
    cpt <- 1
    NnodeLevel <- length(Tree[[level]])
    for(Nnode in 1:NnodeLevel)
    {
      if(Tree[[level]][[Nnode]] == as.character(cptleaf))
      {
        cptleaf <- cptleaf + 1
        next
      }
      L_child <- Tree[[level + 1]][[cpt]] 
      R_child <- Tree[[level + 1]][[cpt + 1]] 
      cpt <- cpt + 2
      adj_list <- rbind(
        adj_list, 
        cbind(
          Tree[[level]][[Nnode]], 
          c(L_child, R_child),
          c("-","+")
        )
      )
    }
  }
  g <- graph.data.frame(data.frame(parent=as.character(adj_list[,1]), 
                                   node=as.character(adj_list[,2]),
                                   text=adj_list[,3])
  )
  E(g)$label.cex <- Ecex
  E(g)$color <- Ecolor
  V(g)$label.cex <- Vcex
  V(g)$color <- Vcolor
  igraph::plot.igraph(g, layout = igraph::layout_as_tree(g), 
                      edge.label=E(g)$text, ...)
  
}
