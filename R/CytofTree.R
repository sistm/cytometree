#' Binary tree algorithm for cytometry data analysis.
#'
#'@param M A matrix of size n x p containing cytometry measures
#'of n cells on p markers.
#'
#'@param minleaf An integer indicating the minimum number of cells
#'per population. Default is \code{1}.
#'
#'@param t A real positive-or-null number used for comparison with
#'the normalized AIC computed at each node of the tree.
#'A higher value limits the height of the tree.
#'
#'@param verbose A logical controlling if a text progress bar is displayed 
#'during the execution of the algorithm. By default is TRUE.
#'
#'@param force_first_markers a vector of index to split the data on first.
#'This argument is used in the semi-supervised setting, forcing the algorithm to consider 
#'those markers first, in the order they appear in this \code{force_first_markers} vector, 
#'and forcing the split at every node. Default is \code{NULL}, in which case
#'the clustering algorithm is unsupervised.
#'
#'@return An object of class 'cytomeTree' providing a partitioning
#'of the set of n cells.
#'\itemize{
#'\item{\code{annotation}}{ A \code{data.frame} containing the annotation of each 
#'cell population underlying the tree pattern.}
#'\item{\code{labels}}{ The partitioning of the set of n cells.}
#'\item{\code{M}}{ The input matrix.}
#'\item{\code{mark_tree}}{ A two level list containing markers used
#'for node splitting.}
#' }
#'
#'@details The algorithm is based on the construction of a binary tree,
#'the nodes of which are subpopulations of cells. At each node,
#'observed cells and markers are modeled by both a family of normal
#'distributions and a family of bi-modal normal mixture distributions.
#'Splitting is done according to a normalized difference of AIC between
#'the two families.
#'@author Chariff Alkhassim, Boris Hejblum
#'
#'@export
#'
#'@examples
#'head(DLBCL)

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
    if(any(!(num_col%in%1:ncol(M)))){
      stop(paste0("'num_col' should contain columns indices from 1 to ", ncol(M)))
    }else{
      M <- AsinhTransformation(M, num_col)
    }
  }
  
  BT <- Cytof_BinaryTree(M, floor(minleaf), t, verbose, force_first_markers)
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
