#'E-M algorithm for cytoftree.
#'
#'@importFrom mclust Mclust
#'
#'@importFrom stats var
#'
#'@importFrom methods is
#'
#'@keywords internal
CytofEM <- function(M, indices, minleaf, level, t, force_marker = NULL){
  
  if(!methods::is(M, "matrix")){
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
  
  aic_norm_old <- -Inf
  parameters  <- list()
  mark_not_dis <- c()
  child <- list()
  
  if(is.null(force_marker)){
    for(j in 1:p){
      flag_uni <- 0
      M_j <- M[,j]
      
      ind_zero <- which(M_j==0) # A
      ind_no_zero <- which(M_j!=0) # A
      
      n_j <- length(ind_no_zero) # A
      if(!var(M_j) | is.na(var(M_j[ind_no_zero]))){ # C
        mark_not_dis <- append(mark_not_dis, colnames(M)[j])
        next() 
      }
      mc_uni <- Mclust(M_j[ind_no_zero], G=1, verbose = FALSE) # C
      mc_mix <- Mclust(M_j[ind_no_zero], G=2, modelNames = "E", verbose = FALSE) # C
      
      ind1_null <- M_j[ind_zero] == 0 # C
      
      ind1 <- mc_mix$classification == 1
      ind2 <- mc_mix$classification == 2
      
      if(length(which(ind1)) < minleaf | length(which(ind2)) < minleaf | is.null(mc_mix)){
        #if(length(ind1) < minleaf | length(ind2) < minleaf | is.null(mc_mix)){
        mark_not_dis <- append(mark_not_dis, colnames(M)[j])
        next()
      }
      
      aic_uni <- 2*mc_uni$df - 2*mc_uni$loglik
      aic_mix <- 2*mc_mix$df - 2*mc_mix$loglik
      aic_norm_new <- (aic_uni - aic_mix)/n_j # C
      
      if(flag_uni | aic_norm_new < t){
        mark_not_dis <- append(mark_not_dis, colnames(M)[j])
      }else{
        
        M1 <- M_j[names(which(ind1==T))] # C
        M2 <- M_j[names(which(ind2==T))] # C
        
        if (mean(M1) > mean(M2)){ # A
          ind1_all <- c(ind1_null, ind2) # A
        }else{ # A
          ind1_all <- c(ind1_null, ind1) # A
        } # A
        
        ind1_all <- ind1_all[order(as.numeric(names(ind1_all)))] # A
        ind2_all <- !ind1_all # A
        
        label <- ind1_all # C
        label[ind1_all] <- 0 # A
        label[ind2_all] <- 1 # A
        
        mean_M1 <- mean(M1)
        mean_M2 <- mean(M2)
        pi_M1 <- length(M1)/n_j # C
        pi_M2 <- 1 - pi_M1
        if(mean_M1 > mean_M2){
          
          temparameters <- cbind.data.frame("aic_norm" = aic_norm_new, 
                                            "marker" = colnames(M)[j], 
                                            "mean_M1" = mean_M2, 
                                            "mean_M2" = mean_M1,
                                            "Var1" = stats::var(M2), 
                                            "Var2" = stats::var(M1), 
                                            "pi1" = pi_M2, 
                                            "pi2" = pi_M1)
        }else{
          
          temparameters <- cbind.data.frame("aic_norm" = aic_norm_new, 
                                            "marker" = colnames(M)[j], 
                                            "mean_M1" = mean_M1, 
                                            "mean_M2" = mean_M2,
                                            "Var1" = stats::var(M1), 
                                            "Var2" = stats::var(M2), 
                                            "pi1" = pi_M1, 
                                            "pi2" = pi_M2)
        }
        if(aic_norm_old < aic_norm_new){
          child$L <-  indices[label == 0]
          child$R <-  indices[label == 1]
          aic_norm_old <- aic_norm_new
        }
        parameters <- rbind.data.frame(parameters, temparameters)
      }
    }
  }else{
    #split on force_marker
    force_marker_index <- which(colnames(M) == force_marker)
    mark_not_dis <- colnames(M)[-force_marker_index]
    
    M_j <- M[, force_marker_index]
    
    ind_zero <- which(M_j==0) # A
    ind_no_zero <- which(M_j!=0) # A
    n_j <- length(ind_no_zero) # A
    
    if(length(ind_no_zero)<2){
      message(paste0("Unable to force split on ", force_marker, " for some node at level", level))
    }else{
      
      mc_uni <- Mclust(M_j[ind_no_zero], G=1, verbose = FALSE) # C
      mc_mix <- Mclust(M_j[ind_no_zero], G=2, modelNames = "E", verbose = FALSE) # C
      
      ind1_null <- M_j[ind_zero] == 0 # C
      
      ind1 <- mc_mix$classification == 1
      ind2 <- mc_mix$classification == 2
      if(sum(ind1)<1 | sum(ind2)<1){
        message(paste0("Unable to force split on ", force_marker, " for some node at level", level))
      }else{
        
        M1 <- M_j[names(which(ind1==T))] # C
        M2 <- M_j[names(which(ind2==T))] # C
        
        aic_uni <- 2*mc_uni$df - 2*mc_uni$loglik
        aic_mix <- 2*mc_mix$df - 2*mc_mix$loglik
        aic_norm_new <- (aic_uni - aic_mix)/n_j # C
        
        if (mean(M1) > mean(M2)){ # A
          ind1_all <- c(ind1_null, ind2) # A
        }else{ # A
          ind1_all <- c(ind1_null, ind1) # A
        } # A
        
        ind1_all <- ind1_all[order(as.numeric(names(ind1_all)))] # A
        ind2_all <- !ind1_all # A
        
        label <- ind1_all # C
        label[ind1_all] <- 0 # A
        label[ind2_all] <- 1 # A
        
        mean_M1 <- mean(M1)
        mean_M2 <- mean(M2)
        pi_M1 <- length(M1)/n_j # C
        pi_M2 <- 1 - pi_M1
        if(mean_M1 > mean_M2){
          
          temparameters <- cbind.data.frame("aic_norm" = aic_norm_new, 
                                            "marker" = force_marker, 
                                            "mean_M1" = mean_M2, 
                                            "mean_M2" = mean_M1,
                                            "Var1" = stats::var(M2), 
                                            "Var2" = stats::var(M1), 
                                            "pi1" = pi_M2, 
                                            "pi2" = pi_M1)
        }else{
          
          temparameters <- cbind.data.frame("aic_norm" = aic_norm_new, 
                                            "marker" = force_marker, 
                                            "mean_M1" = mean_M1, 
                                            "mean_M2" = mean_M2,
                                            "Var1" = stats::var(M1), 
                                            "Var2" = stats::var(M2), 
                                            "pi1" = pi_M1, 
                                            "pi2" = pi_M2)
        }
        
        child$L <-  indices[label == 0]
        child$R <-  indices[label == 1]
        parameters <- rbind.data.frame(parameters, temparameters)
      }
    }
  }
  
  nnrowpara <- nrow(parameters)
  
  ans <- NULL
  if(is.null(nnrowpara)){
    ans <- list("mark_not_dis" = mark_not_dis)
  }else{
    if(nnrowpara > 1){
      parameters <- parameters[order(parameters[, "aic_norm"], decreasing = TRUE), ]
    }
    ans <- list("mark_not_dis" = mark_not_dis, "child" = child,
                "nAIC" = parameters[, 1], "ind" = as.character(parameters[, 2]),
                "mu1"= parameters[1, 3], "mu2"= parameters[1, 4],
                "Var1" = parameters[1, 5], "Var2" = parameters[1, 6],
                "pi1"= parameters[1, 7], "pi2" = parameters[1, 8])
  }
  
  return(ans)
}