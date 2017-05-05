#' TODO
#' 
#'@param 
#'
#'@author Chariff Alkhassim
#'
#'@export 
# 


Annotation<- function(CytomeTreeObj, K3markers, plot)
{
  if(class(CytomeTreeObj) != "CytomeTree")
  {
    stop("CytomeTreeObj must be of class CytomeTree")
  }
  if(!is.null(K3markers))
  {
    if(class(K3markers)!="character")
    {
      stop("K3markers must be of class character")
    }
  } 
  M <- CytomeTreeObj$M
  labels <- CytomeTreeObj$labels
  lc <- LeafCenters(CytomeTreeObj)
  len_lab <- length(labels)
  dlc <- dim(lc)
  n <- dlc[1]
  p <- dlc[2]
  leaves <- lc[,p]
  cnames <- colnames(lc[,1:(p-1)])
  combinations <- cbind(matrix(0, ncol = (p-1), nrow = n), 1:n)
  if(n == 1) 
  {
    stop("CytomeTree didn't find any population")
  }
  else 
  {
    for(j in 1:(p-1))
    {
      ExpressLevels <- 2
      if(any(cnames[j] == K3markers))
      {
        ExpressLevels <- 3
      }
      leavesSort <- leaves[sort(lc[,j], index.return = TRUE)$ix]
      M_j <- M[,j]
      leavesSort_ <- leavesSort
      partitions2gr <- Partition2gr(n)
      Kmeans2 <- KmeansOPT(partitions2gr, leavesSort, labels, M_j, K = 2)
      partwin2gr <- partitions2gr[[Kmeans2$ind]]
      tempclass_neg.2  <- leavesSort[which(partwin2gr == 1)]
      tempclass_pos.2  <- leavesSort[which(partwin2gr == 2)]
      tind1.2 <- which(labels%in%tempclass_neg.2)
      tind2.2 <- which(labels%in%tempclass_pos.2)
      partitions3gr <- Partition3gr(n)
      Kmeans3 <- KmeansOPT(partitions3gr, leavesSort, labels, M_j, K = 3)
      partwin3gr <- partitions3gr[[which.max(Kmeans3$ind)]]
      tempclass_neg.3  <- leavesSort[which(partwin3gr == 1)]
      tempclass_pos.3  <- leavesSort[which(partwin3gr == 2)]
      tempclass_dpos.3 <- leavesSort[which(partwin3gr == 3)]
      tind1.3 <- which(labels%in%tempclass_neg.3)
      tind2.3 <- which(labels%in%tempclass_pos.3)
      tind3.3 <- which(labels%in%tempclass_dpos.3)
      if(ExpressLevels == 2)
      {   
        combinations[tempclass_pos.2, j] <- 1
        if(plot)
        {
          Expression <- rep(1, len_lab)
          Expression[labels%in%tempclass_neg.2] <- 2      
          dfbox <- data.frame(Leaves = factor(labels, levels = 
                                                as.character(leavesSort)), 
                              Fluorescence = M[,j], 
                              Expression = as.factor(Expression))
          p <- ggplot(dfbox, aes(Leaves, 
                                 Fluorescence, 
                                 fill = Expression ))
          suppressWarnings(print(p + ggtitle(cnames[j]) + 
                                   geom_boxplot(outlier.shape = NA, alpha = 1/3)+
                                   scale_fill_manual(values = c("red","blue"),
                                                     name = "Annotation",
                                                     labels = c("Hi","Low"))))
        }
      }
      else if(ExpressLevels == 3)
      {
        combinations[tempclass_pos.3, j] <- 1
        combinations[tempclass_dpos.3, j] <- 2
        if(plot)
        {        
          Expression <- rep(1, len_lab)
          Expression[labels%in%tempclass_pos.3] <- 2
          Expression[labels%in%tempclass_dpos.3] <- 3
          
          
          dfbox <- data.frame(Leaves = factor(labels, levels = 
                                                as.character(leavesSort)), 
                              Fluorescence = M[,j], 
                              Expression = as.factor(Expression))
          p <- ggplot(dfbox, aes(Leaves, 
                                 Fluorescence, 
                                 fill = Expression ))
          suppressWarnings(print(p + ggtitle(cnames[j]) + 
                                   geom_boxplot(outlier.shape = NA, 
                                                alpha = 1/3)+
                                   scale_fill_manual(values = 
                                                       c("blue","green","red"),
                                                     name = "Annotation",
                                                     labels = 
                                                       c("Low","Hi","Hi+"))))
        }
      }
    } 
  }
  
  tblabels <- table(labels)
  combinations <- cbind(combinations, table(labels), round(tblabels/len_lab,4))
  colnames(combinations) <- c(cnames, "leaves", "count", "prop")
  as.data.frame(combinations)
}




