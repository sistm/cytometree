#' Gridsearch
#' 
#' @export

CytGridsearch <- function(M, indices, minleaf, level, t, force_marker = NULL, cytof = FALSE){
  
  if(class(M)!="matrix"){
    warning("M is not a matrix, this should not happen\n Please let the package maintainer know that ', drop=FALSE' is probably missing somewhere...")
    M <- as.matrix(M)
  }
  
  n <- nrow(M)
  p <- ncol(M)
  
  if(n <= minleaf ){
    return(list("mark_not_dis" = colnames(M)))
  }
  
  nEMdegenerate <- 5
  
  if(minleaf < nEMdegenerate){
    minleaf <- nEMdegenerate
  }
  
  parameters <- c()
  child <- list()
  aic_n_old <- -Inf
  
  mark_not_dis <- c()
  mark_split <- c()
  
  #browser()
  
  if(is.null(force_marker)){  
    if (cytof){
      for (j in 1:p){
        
        M_j <- M[,j]
        
        ind_zero <- which(M_j==0)
        ind_no_zero <- which(M_j!=0)
        
        if(!var(M_j) | is.na(var(M_j[ind_no_zero]))){
          mark_not_dis <- append(mark_not_dis, colnames(M)[j])
          next() 
        }
        
        n_j <- length(ind_no_zero)
        
        res1_MLE <- aic_1_gauss(x = M_j[ind_no_zero])
        res2_MLE <- Gridsearch(x = M_j[ind_no_zero], resu_aic_1_gauss = res1_MLE, 
                               iter_max = 15, ntry = 30, mixture = 2)
        
        aic_n <- (res1_MLE$AIC - res2_MLE$AIC)/n_j
        
        if(aic_n < t){
          mark_not_dis <- append(mark_not_dis, colnames(M)[j])
          next()
        }else{
          
          ind1_null <- M_j[ind_zero] == 0
          
          clust1 <- res2_MLE$cluster == 1
          clust2 <- res2_MLE$cluster == 2
          
          if (res2_MLE$mu1 > res2_MLE$mu2){
            ind1_all <- c(ind1_null, clust2)
          }else{
            ind1_all <- c(ind1_null, clust1)
          }
          
          ind1_all <- ind1_all[order(as.numeric(names(ind1_all)))]
          ind2_all <- !ind1_all
          
          label <- ind1_all
          label[ind1_all] <- 0 
          label[ind2_all] <- 1
          
          if(length(which(ind1_all)) < minleaf | length(which(ind2_all)) < minleaf){
            mark_not_dis <- append(mark_not_dis, colnames(M)[j])
            next()
          }
          
          if (res2_MLE$mu1 > res2_MLE$mu2){
            
            temp_parameters <- cbind.data.frame("marker" = as.character(colnames(M)[j]),
                                                "aic_n" = aic_n,
                                                "mu1" = res2_MLE$mu2,
                                                "mu2" = res2_MLE$mu1,
                                                "sigma1" = res2_MLE$s2,
                                                "sigma2" = res2_MLE$s1,
                                                "p" = 1 - res2_MLE$p)
            
          }else{
            
            temp_parameters <- cbind.data.frame("marker" = as.character(colnames(M)[j]),
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
      
    }else{
      
      for (j in 1:p){
        
        M_j <- M[,j]
        
        if(!var(M_j)){
          mark_not_dis <- append(mark_not_dis, colnames(M)[j])
          next() 
        }
        
        n_j <- length(M_j)
        
        res1_MLE <- aic_1_gauss(x = M_j)
        res2_MLE <- Gridsearch(x = M_j, resu_aic_1_gauss = res1_MLE, 
                               iter_max = 15, ntry = 30, mixture = 2)
        
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
            
            temp_parameters <- cbind.data.frame("marker" = as.character(colnames(M)[j]),
                                                "aic_n" = aic_n,
                                                "mu1" = res2_MLE$mu2,
                                                "mu2" = res2_MLE$mu1,
                                                "sigma1" = res2_MLE$s2,
                                                "sigma2" = res2_MLE$s1,
                                                "p" = 1 - res2_MLE$p)
            
          }else{
            
            label[clust1] <- 0
            label[clust2] <- 1
            
            temp_parameters <- cbind.data.frame("marker" = as.character(colnames(M)[j]),
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
    }
    
  }else{
    #split on force_marker
    force_marker_index <- which(colnames(M) == force_marker)
    mark_not_dis <- colnames(M)[-force_marker_index]
    
    if (cytof){
      
      M_j <- M[,force_marker_index]
      
      ind_zero <- which(M_j==0)
      ind_no_zero <- which(M_j!=0)
      
      if(!var(M_j) | is.na(var(M_j[ind_no_zero]))){
        message(paste0("Unable to force split on ", force_marker, " for some node at level", level))
      }else{
        
        n_j <- length(ind_no_zero)
        
        res1_MLE <- aic_1_gauss(x = M_j[ind_no_zero])
        res2_MLE <- Gridsearch(x = M_j[ind_no_zero], resu_aic_1_gauss = res1_MLE, 
                               iter_max = 15, ntry = 30, mixture = 2)
        
        aic_n <- (res1_MLE$AIC - res2_MLE$AIC)/n_j
        
        ind1_null <- M_j[ind_zero] == 0
        
        clust1 <- res2_MLE$cluster == 1
        clust2 <- res2_MLE$cluster == 2
        
        if (res2_MLE$mu1 > res2_MLE$mu2){
          ind1_all <- c(ind1_null, clust2)
        }else{
          ind1_all <- c(ind1_null, clust1)
        }
        
        ind1_all <- ind1_all[order(as.numeric(names(ind1_all)))]
        ind2_all <- !ind1_all
        
        if(sum(ind1_all) < 1 | sum(ind2_all) < 1){
          message(paste0("Unable to force split on ", force_marker, " for some node at level", level))
        }else{
          
          label <- ind1_all
          label[ind1_all] <- 0 
          label[ind2_all] <- 1
          
          if (res2_MLE$mu1 > res2_MLE$mu2){
            
            temp_parameters <- cbind.data.frame("marker" = force_marker,
                                                "aic_n" = aic_n,
                                                "mu1" = res2_MLE$mu2,
                                                "mu2" = res2_MLE$mu1,
                                                "sigma1" = res2_MLE$s2,
                                                "sigma2" = res2_MLE$s1,
                                                "p" = 1 - res2_MLE$p)
            
          }else{
            
            temp_parameters <- cbind.data.frame("marker" = force_marker,
                                                "aic_n" = aic_n,
                                                "mu1" = res2_MLE$mu1,
                                                "mu2" = res2_MLE$mu2,
                                                "sigma1" = res2_MLE$s1,
                                                "sigma2" = res2_MLE$s2,
                                                "p" = res2_MLE$p)
            
          }
          
          parameters <- rbind(parameters, temp_parameters)
          
          child$L <- indices[label == 0]
          child$R <- indices[label == 1]
          
          mark_split <- force_marker
          
        }
      }
      
    }else{
      
      M_j <- M[,force_marker_index]
      
      n_j <- length(M_j)
      
      res1_MLE <- aic_1_gauss(x = M_j)
      res2_MLE <- Gridsearch(x = M_j, resu_aic_1_gauss = res1_MLE, 
                             iter_max = 15, ntry = 30, mixture = 2)
      
      clust1 <- res2_MLE$cluster == 1
      clust2 <- res2_MLE$cluster == 2
      
      if(sum(clust1)<1 | sum(clust2)<1){
        message(paste0("Unable to force split on ", force_marker, " for some node at level", level))
      }else{
        
        aic_n <- (res1_MLE$AIC - res2_MLE$AIC)/n_j
        
        label <- res2_MLE$cluster
        
        if (res2_MLE$mu1 > res2_MLE$mu2){
          
          label[clust1] <- 1
          label[clust2] <- 0
          
          temp_parameters <- cbind.data.frame("marker" = force_marker,
                                              "aic_n" = aic_n,
                                              "mu1" = res2_MLE$mu2,
                                              "mu2" = res2_MLE$mu1,
                                              "sigma1" = res2_MLE$s2,
                                              "sigma2" = res2_MLE$s1,
                                              "p" = 1 - res2_MLE$p)
          
        }else{
          
          label[clust1] <- 0
          label[clust2] <- 1
          
          temp_parameters <- cbind.data.frame("marker" = force_marker,
                                              "aic_n" = aic_n,
                                              "mu1" = res2_MLE$mu1,
                                              "mu2" = res2_MLE$mu2,
                                              "sigma1" = res2_MLE$s1,
                                              "sigma2" = res2_MLE$s2,
                                              "p" = res2_MLE$p)
          
        }
        
        parameters <- rbind(parameters, temp_parameters)
        
        child$L <- indices[label == 0]
        child$R <- indices[label == 1]
        mark_split <- force_marker
        
      }
    }
  }  
  
  if (is.null(mark_split)){
    
    return(list("mark_not_dis" = mark_not_dis))
    
  }else{
    
    if (nrow(parameters) > 1){
      parameters[order(parameters[, "aic_n"], decreasing = TRUE), ]
    }
    
    return(list("mark_not_dis" = mark_not_dis, "child" = child,
                "nAIC" = parameters$aic_n, "ind" = as.character(parameters$marker),
                "mu1"= parameters[1,"mu1"], "mu2"= parameters[1,"mu2"],
                "Var1" = parameters[1,"sigma1"], "Var2" = parameters[1,"sigma2"],
                "pi1"= parameters[1,"p"], "pi2" = 1 - parameters[1,"p"]))
    
  }
  
}
