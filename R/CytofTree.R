CytofTree <- function(M, minleaf = 1, t = .1, verbose = TRUE, 
                      force_first_markers = NULL, transformation = c("asinh", "none"), 
                      num_col = NULL){
  if(class(M) == "data.frame"){
    M <- as.matrix(M)
  }
  if(class(M) != "matrix" | mode(M) != "numeric"){
    stop("M should be a numeric matrix")
  }
  n <- nrow(M)
  if(minleaf >= n)
  {
    stop("minleaf is superior to n.")
  }
  p <- ncol(M)
  if(p > n){
    stop("p is superior to n.")
  }
  if(any(is.na(M)))
  {
    stop("M contains NAs.")
  }
  
  if(is.null(colnames(M))){
    colnames(M) <- paste0(rep("M",p), 1:p)
  }
  
  if(!is.null(force_first_markers)){
    if(any(!(force_first_markers %in% colnames(M)))){
      cat("'force_first_markers' are not all in M colnames:")
      colnames(M)
      stop()
    }
  }
  
  if(any(!transformation%in%c("asinh", "none"))){
    stop("The only allowed values for transformation are 'asinh' or 'none'")
  }
  
  if(any(transformation%in%c("asinh"))){
    if(any(!(1:ncol(M)%in%num_col))){
      stop(paste0("'num_col' should contain columns indices from 1 to ", ncol(M)))
    }else{
      M <- AsinhTransformation(M, num_col)
    }
  }
  
  BT <- cytof_BinaryTree(M, floor(minleaf), t, verbose, force_first_markers)
  annotation <- TreeAnnot(BT$labels, BT$combinations)
  Tree <- list("M" = M, "labels" = BT$labels,
               "pl_list"= BT$pl_list, "t"= t,
               "mark_tree" = BT$mark_tree,
               "annotation" = annotation,
               "transformation" = transformation,
               "num_col" = num_col)
  class(Tree) <- "CytomeTree"
  cat("It works !", "\n")
  return(Tree)
}
