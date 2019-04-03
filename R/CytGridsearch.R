#' Gridsearch
#' 
#' @export

CytGridsearch <- function(M, indices, minleaf, level, t, force_marker = NULL, cytof = FALSE){
  
  n <- nrow(M)
  p <- ncol(M)
  parameters <- c()
  child <- list()
  aic_n_old <- -Inf
  
  mark_not_dis <- c()
  mark_split <- c()
  
  #browser()
  
  if (cytof){
    
  }else{
    for (j in 1:p){
      
      M_j <- M[,j]
      
      n_j <- length(M_j)
      
      res1_MLE <- aic_1_gauss(M_j)
      res2_MLE <- Gridsearch(x = M_j, resu_aic_1_gauss = res1_MLE, iter_max = 100, ntry = 30)
      
      aic_n <- (res1_MLE$AIC - res2_MLE$AIC)/n_j
      
      if(aic_n < t){
        mark_not_dis <- append(mark_not_dis, colnames(M)[j])
        next()
      }else{
        
        label <- res2_MLE$cluster
        clust1 <- res2_MLE$cluster == 1
        clust2 <- res2_MLE$cluster == 2
        
        if(length(which(clust1)) < minleaf | length(which(clust2)) < minleaf){
          mark_not_dis <- append(mark_not_dis, colnames(M)[j])
          next()
        }
        
        if (res2_MLE$mu1 > res2_MLE$mu2){
          
          label[clust1] <- 1
          label[clust2] <- 0
          
          temp_parameters <- cbind.data.frame("marker" = colnames(M)[j],
                                              "aic_n" = aic_n,
                                              "mu1" = res2_MLE$mu2,
                                              "mu2" = res2_MLE$mu1,
                                              "sigma1" = res2_MLE$s2,
                                              "sigma2" = res2_MLE$s1,
                                              "p" = 1 - res2_MLE$p)
          
        }else{
          
          label[clust1] <- 0
          label[clust2] <- 1
          
          temp_parameters <- cbind.data.frame("marker" = colnames(M)[j],
                                              "aic_n" = aic_n,
                                              "mu1" = res2_MLE$mu1,
                                              "mu2" = res2_MLE$mu2,
                                              "sigma1" = res2_MLE$s1,
                                              "sigma2" = res2_MLE$s2,
                                              "p" = res2_MLE$p)
          
        }
        
        parameters <- rbind(parameters, temp_parameters)
        
        if (aic_n > aic_n_old){
          
          aic_n_old <- aic_n
          
          child$L <- indices[label == 0]
          child$R <- indices[label == 1]
          
          mark_split <- as.character(colnames(M)[j])
        }
      }
    }
    
    if (is.null(mark_split)){
      
      return(list("mark_not_dis" = mark_not_dis))
      
    }else{
      
      best_parameter <- parameters[which.max(parameters[,"aic_n"]),]
      
      return(list("mark_not_dis" = mark_not_dis, "child" = child,
                  "nAIC" = best_parameter$aic_n, "ind" = as.character(best_parameter$marker),
                  "mu1"= best_parameter$mu1, "mu2"= best_parameter$mu2,
                  "Var1" = best_parameter$sigma1, "Var2" = best_parameter$sigma2,
                  "pi1"= best_parameter$p, "pi2" = 1 - best_parameter$p))
    }
  }
  
  
  
}
