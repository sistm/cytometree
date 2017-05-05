#' TODO
#'@param 
#'
#'@author Chariff Alkhassim
#'
#'@export 
# 
plot_graph <- function(CytomeTreeObj, Ecex = 1, Ecolor = 8,
                       Vcex = .8, Vcolor = 0, ...)
{
  if(class(CytomeTreeObj) != "CytomeTree")
  {
    stop("CytomeTreeObj must be of class CytomeTree")
  }
  Signtree <- unlist(CytomeTreeObj$Signtree)
  Tree <- CytomeTreeObj$mark_tree
  Tree_level <- length(Tree)
  adj_list <- c()
  for(level in 1:(Tree_level - 1))
  {
    cpt <- 1
    NnodeLevel <- length(Tree[[level]])
    for(Nnode in 1:NnodeLevel)
    {
      L_son <- Tree[[level + 1]][[cpt]] 
      R_son <- Tree[[level + 1]][[cpt + 1]] 
      cpt <- cpt + 2
      adj_list <- rbind(adj_list, cbind(Tree[[level]][[Nnode]], 
                                        c(L_son, R_son)))
    }
  }
  rm <- which(rowSums(is.na(adj_list))>0)
  if(length(rm))
  {
    adj_list_ <- adj_list[-c(rm),] 
  }
  else
  {
    adj_list_ <- adj_list
  }
  rm2 <- which(is.na(Signtree))
  if(length(rm2))
  {
    Signtree_ <- Signtree[-c(rm2)]
  }
  else
  {
    Signtree_ <- Signtree
  }
  g <- graph.data.frame(data.frame(parent=as.character(adj_list_[,1]), 
                                   node=as.character(adj_list_[,2]),
                                   text=Signtree_[-c(1)]))
  E(g)$label.cex <- Ecex
  E(g)$color <- Ecolor
  V(g)$label.cex <- Vcex
  V(g)$color <- Vcolor  
  plot(g, layout = layout.reingold.tilford(g),edge.label=E(g)$text, ...)
}
