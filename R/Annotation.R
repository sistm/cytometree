#' Annotates cell populations found using CytomeTree. 
#' 
#'@param CytomeTreeObj An object of class CytomeTree.
#'
#'@param K2markers A vector of class character where the names of 
#'the markers for which 2 levels of expression are sought can be specified.
#'Default is \code{NULL} i.e. unsupervised. 
#'
#'@param K3markers A vector of class character where the names of 
#'the markers for which 3 levels of expression are sought can be specified.
#'Default is \code{NULL} i.e. unsupervised. 
#'
#'@param plot A logical value indicating whether or not to plot the 
#'partitioning in 1, 2 or 3 groups for each marker. Default is \code{TRUE}.
#'
#'@param t A real positive-or-null number used for comparison with
#'the normalized AIC computed to compare the fits of the marginal distributions
#'obtained by one normal distribution and by a mixture of two or three normal.
#'For markers used in the tree, the algorithm compares the
#'fits obtained by a mixture of two and three normal distributions.
#'Default value is .2. A higher value leads to a smaller number of
#'expression levels per marker.
#'
#'@return A \code{data.frame} containing the annotation of each 
#'cell population. 
#'
#'@details The algorithm is set to find the partitioning in 1, 2 or
#'3 groups of cell populations found using CytomeTree. In an unsupervised mode,
#'it minimizes the within-leaves sum of squares of the observed values on each
#'marker and computes the normalized AIC to compare the fits of the marginal 
#'distributions obtained by one normal distribution and by a mixture of two
#'or three normal.For markers used in the tree, the algorithm compares the
#'fits obtained by a mixture of two and three normal distributions.
#'
#'@author Chariff Alkhassim
#'
#'@import ggplot2 graphics mclust
#'
#'@importFrom stats dnorm sd
#'
#'@export
Annotation <- function(CytomeTreeObj, K2markers = NULL, 
                       K3markers = NULL, plot = TRUE, t = 0.2)
{
  if(class(CytomeTreeObj) != "CytomeTree")
  {
    stop("CytometreeObj must be of class CytomeTree.")
  }
  if(!is.null(K3markers))
  {
    if(class(K3markers)!="character")
    {
      stop("K3markers must be of class character.")
    }
  }
  if(is.null(t))
  {
    t <- CytomeTreeObj$t
  }
  M <- CytomeTreeObj$M
  labels <- CytomeTreeObj$labels
  lc <- LeavesCenters(CytomeTreeObj)
  len_lab <- length(labels)
  dlc <- dim(lc)
  n <- dlc[1]
  p <- dlc[2]
  leaves <- lc[,p]
  cnames <- colnames(lc[,1:(p-1)])  
  combinations <- cbind(matrix(0, ncol = (p-1), nrow = n), 1:n)
  if(n == 1) 
  {
    return(cat("CytomeTree found a single population.\n"))
  }
  else 
  {
    indToAuto <- colSums(apply(CytomeTreeObj$annotation[,1:(p-1)],2,is.na)) == n
    FlagMArkerinTree <- seq_len((p-1))[as.logical(abs(indToAuto*1-1))]
    for(j in 1:(p-1))
    {
      leavesSort <- leaves[sort(lc[,j], index.return = TRUE)$ix]
      M_j <- M[,j]
      leavesSort_ <- leavesSort
      if(any(cnames[j] == K3markers))
      {
        partitions3gr <- Partition3gr(n)
        Kmeans3 <- KmeansOPT(partitions3gr, leavesSort, labels, M_j, K = 3)
        partwin3gr <- partitions3gr[[Kmeans3$ind]]
        tempclass_neg.3  <- leavesSort[partwin3gr == 1]
        tempclass_pos.3  <- leavesSort[partwin3gr == 2]
        tempclass_dpos.3 <- leavesSort[partwin3gr == 3]
        tind1.3 <- labels%in%tempclass_neg.3
        tind2.3 <- labels%in%tempclass_pos.3
        tind3.3 <- labels%in%tempclass_dpos.3
        combinations[tempclass_pos.3, j] <- 1
        combinations[tempclass_dpos.3, j] <- 2
        if(plot)
        {   
          Expression <- rep(1, len_lab)
          Expression[tind2.3] <- 2
          Expression[tind3.3] <- 3
          dfbox <- data.frame(Leaves = factor(labels, levels = 
                                                as.character(leavesSort)), 
                              Fluorescence = M[,j], 
                              Expression = as.factor(Expression))
          tempdf <-factor(dfbox[,3], levels=c(3,2,1))
          dfbox[,3] <- tempdf
          p <- ggplot2::ggplot(dfbox, ggplot2::aes_string("Leaves", 
                                                          "Fluorescence", 
                                                          fill="Expression"))+
            ggplot2::xlab("Populations")+
            ggplot2::theme(axis.title=element_text(size=15),
                           axis.text=element_text(size=12,face = "bold"),
                           legend.title=element_text(face="bold"))
          
          suppressWarnings(print(p + ggplot2::ggtitle(cnames[j]) + 
                                   ggplot2::geom_boxplot(outlier.shape=NA, 
                                                         alpha=1)+
                                   ggplot2::scale_fill_manual(values=c("tomato1",
                                                                       "green",
                                                                       "cyan3"),
                                                              name="Annotation",
                                                              labels=c("++",
                                                                       "+",
                                                                       "-")
                                   )
          ))
        }
      }
      else if(any(cnames[j] == K2markers))
      {
        partitions2gr <- Partition2gr(n)
        Kmeans2 <- KmeansOPT(partitions2gr, leavesSort, labels, M_j, K = 2)
        partwin2gr <- partitions2gr[[Kmeans2$ind]]
        tempclass_neg.2  <- leavesSort[partwin2gr == 1]
        tempclass_pos.2  <- leavesSort[partwin2gr == 2]
        tind1.2 <- labels%in%tempclass_neg.2
        tind2.2 <- labels%in%tempclass_pos.2
        combinations[tempclass_pos.2, j] <- 1
        if(plot)
        {
          Expression <- rep(1, len_lab)
          Expression[tind1.2] <- 2     
          
          dfbox <- data.frame(Leaves = factor(labels, levels = 
                                                as.character(leavesSort)), 
                              Fluorescence = M[,j], 
                              Expression = as.factor(Expression))
          p <- ggplot2::ggplot(dfbox, ggplot2::aes_string("Leaves", 
                                                          "Fluorescence", 
                                                          fill = "Expression"))+
            ggplot2::xlab("Populations")+
            ggplot2::theme(axis.title=element_text(size=15),
                           axis.text=element_text(size=12,face = "bold"),
                           legend.title=element_text(face="bold"))
          suppressWarnings(print(p + ggplot2::ggtitle(cnames[j]) + 
                                   ggplot2::geom_boxplot(outlier.shape = NA, 
                                                         alpha = 1)+
                                   ggplot2::scale_fill_manual(values=c("tomato1",
                                                                       "cyan3"),
                                                              name="Annotation",
                                                              labels=c("+",
                                                                       "-")
                                   )
          ))
        }
      }
      else
      {
        partitions2gr <- Partition2gr(n)
        Kmeans2 <- KmeansOPT(partitions2gr, leavesSort, labels, M_j, K = 2)
        partwin2gr <- partitions2gr[[Kmeans2$ind]]
        tempclass_neg.2  <- leavesSort[partwin2gr == 1]
        tempclass_pos.2  <- leavesSort[partwin2gr == 2]
        tind1.2 <- labels%in%tempclass_neg.2
        tind2.2 <- labels%in%tempclass_pos.2
        partitions3gr <- Partition3gr(n)
        Kmeans3 <- KmeansOPT(partitions3gr, leavesSort, labels, M_j, K = 3)
        partwin3gr <- partitions3gr[[Kmeans3$ind]]
        tempclass_neg.3  <- leavesSort[partwin3gr == 1]
        tempclass_pos.3  <- leavesSort[partwin3gr == 2]
        tempclass_dpos.3 <- leavesSort[partwin3gr == 3]
        tind1.3 <- labels%in%tempclass_neg.3
        tind2.3 <- labels%in%tempclass_pos.3
        tind3.3 <- labels%in%tempclass_dpos.3
        comp1.2 <- M_j[tind1.2]
        comp2.2 <- M_j[tind2.2]
        comp1.3 <- M_j[tind1.3]
        comp2.3 <- M_j[tind2.3]
        comp3.3 <- M_j[tind3.3]
        mc2 <- Mclust(M_j, G=2, modelNames = "E", verbose = FALSE)
        mc3 <- Mclust(M_j, G=3, modelNames = "E", verbose = FALSE)
        aic2 <- 2*mc2$df - 2*mc2$loglik
        aic3 <- 2*mc3$df - 2*mc3$loglik
        aic_norm_23 <- (aic2 - aic3)/len_lab
      
        if(any(FlagMArkerinTree%in%j))
        {
          if(aic_norm_23 > t)
          {
            combinations[tempclass_pos.3, j] <- 1
            combinations[tempclass_dpos.3, j] <- 2
            if(plot)
            {        
              Expression <- rep(1, len_lab)
              Expression[tind2.3] <- 2
              Expression[tind3.3] <- 3
              dfbox <- data.frame(Leaves = factor(labels, levels = 
                                                    as.character(leavesSort)), 
                                  Fluorescence = M[,j], 
                                  Expression = as.factor(Expression))
              tempdf <-factor(dfbox[,3], levels=c(3,2,1))
              dfbox[,3] <- tempdf
              p <- ggplot2::ggplot(dfbox, ggplot2::aes_string("Leaves", 
                                                              "Fluorescence", 
                                                              fill="Expression"))+
                ggplot2::xlab("Populations")+
                ggplot2::theme(axis.title=element_text(size=15),
                               axis.text=element_text(size=12,face = "bold"),
                               legend.title=element_text(face="bold"))
              
              suppressWarnings(print(p + ggplot2::ggtitle(cnames[j]) + 
                                       ggplot2::geom_boxplot(outlier.shape=NA, 
                                                             alpha=1)+
                                       ggplot2::scale_fill_manual(values=c("tomato1",
                                                                           "green",
                                                                           "cyan3"),
                                                                  name="Annotation",
                                                                  labels=c("++",
                                                                           "+",
                                                                           "-")
                                       )
              ))
            }
          }
          else
          {
            combinations[tempclass_pos.2, j] <- 1
            if(plot)
            {
              Expression <- rep(1, len_lab)
              Expression[tind1.2] <- 2      
              dfbox <- data.frame(Leaves = factor(labels, levels = 
                                                    as.character(leavesSort)), 
                                  Fluorescence = M[,j], 
                                  Expression = as.factor(Expression))
              p <- ggplot2::ggplot(dfbox, ggplot2::aes_string("Leaves", 
                                                              "Fluorescence", 
                                                              fill = "Expression"))+
                ggplot2::xlab("Populations")+
                ggplot2::theme(axis.title=element_text(size=15),
                               axis.text=element_text(size=12,face = "bold"),
                               legend.title=element_text(face="bold"))
              
              suppressWarnings(print(p + ggplot2::ggtitle(cnames[j]) + 
                                       ggplot2::geom_boxplot(outlier.shape = NA, 
                                                             alpha = 1)+
                                       ggplot2::scale_fill_manual(values=c("tomato1",
                                                                           "cyan3"),
                                                                  name="Annotation",
                                                                  labels=c("+",
                                                                           "-")
                                       )
              ))
            }
          }
        }
        else
        {
          aic1 <- 4 - 2*sum(stats::dnorm(M_j, mean(M_j), sd(M_j), log = TRUE))
          aic_norm_12 <- (aic1 - aic2)/len_lab
          if(aic_norm_12 > t)
          {
            if(aic_norm_23 > t)
            {
              combinations[tempclass_pos.3, j] <- 1
              combinations[tempclass_dpos.3, j] <- 2
              if(plot)
              {        
                Expression <- rep(1, len_lab)
                Expression[tind2.3] <- 2
                Expression[tind3.3] <- 3
                dfbox <- data.frame(Leaves = factor(labels, levels = 
                                                      as.character(leavesSort)), 
                                    Fluorescence = M[,j], 
                                    Expression = as.factor(Expression))
                tempdf <-factor(dfbox[,3], levels=c(3,2,1))
                dfbox[,3] <- tempdf
                p <- ggplot2::ggplot(dfbox, ggplot2::aes_string("Leaves", 
                                                                "Fluorescence", 
                                                                fill="Expression"))+
                  ggplot2::xlab("Populations")+
                  ggplot2::theme(axis.title=element_text(size=15),
                                 axis.text=element_text(size=12,face = "bold"),
                                 legend.title=element_text(face="bold"))
                
                suppressWarnings(print(p + ggplot2::ggtitle(cnames[j]) + 
                                         ggplot2::geom_boxplot(outlier.shape=NA, 
                                                               alpha=1)+
                                         ggplot2::scale_fill_manual(values=c("tomato1",
                                                                             "green",
                                                                             "cyan3"),
                                                                    name="Annotation",
                                                                    labels=c("++",
                                                                             "+",
                                                                             "-")
                                         )
                ))
              }
            }
            else
            {
              combinations[tempclass_pos.2, j] <- 1
              if(plot)
              {
                Expression <- rep(1, len_lab)
                Expression[tind1.2] <- 2      
                dfbox <- data.frame(Leaves = factor(labels, levels = 
                                                      as.character(leavesSort)), 
                                    Fluorescence = M[,j], 
                                    Expression = as.factor(Expression))
                p <- ggplot2::ggplot(dfbox, ggplot2::aes_string("Leaves", 
                                                                "Fluorescence", 
                                                                fill = "Expression"))+
                  ggplot2::xlab("Populations")+
                  ggplot2::theme(axis.title=element_text(size=15),
                                 axis.text=element_text(size=12,face = "bold"),
                                 legend.title=element_text(face="bold"))
                
                suppressWarnings(print(p + ggplot2::ggtitle(cnames[j]) + 
                                         ggplot2::geom_boxplot(outlier.shape = NA, 
                                                               alpha = 1)+
                                         ggplot2::scale_fill_manual(values=c("tomato1",
                                                                             "cyan3"),
                                                                    name="Annotation",
                                                                    labels=c("+",
                                                                             "-")
                                         )
                ))
              }
            }
          }
          else
          {
            combinations[, j] <- NA
            if(plot)
            {    
              dfbox <- data.frame(Leaves = factor(labels, levels = 
                                                    as.character(leavesSort)), 
                                  Fluorescence = M[,j])
              p <- ggplot2::ggplot(dfbox, ggplot2::aes_string("Leaves", 
                                                              "Fluorescence"))+
                ggplot2::xlab("Populations")+
                ggplot2::theme(axis.title=element_text(size=15),
                               axis.text=element_text(size=12,face = "bold"),
                               legend.title=element_text(face="bold"))
              suppressWarnings(print(p + ggplot2::ggtitle(cnames[j]) + 
                                       ggplot2::geom_boxplot(outlier.shape = NA, 
                                                             alpha = 1)                       
              ))
            } 
          }
        }
      }
    } 
  }
  tblabels <- table(labels)
  combinations <- cbind(combinations, table(labels), round(tblabels/len_lab,4))
  colnames(combinations) <- c(cnames, "leaves", "count", "prop")
  outCombinations <- as.data.frame(combinations[sort(combinations[,"count"],
                                                     TRUE,
                                                     index.return=TRUE)$ix,])  
  out <- list("combinations" = outCombinations, "labels" = labels)
  class(out) <- "Annotation"
  out
}




