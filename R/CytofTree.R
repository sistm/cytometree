#' Binary tree algorithm for mass cytometry data analysis.
#'
#'@param M A matrix of size n x p containing mass cytometry measures
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
#'@param transformation A string indicating the transformation used among \code{asinh} 
#'\code{biexp}, \code{log10} and \code{none}. Default is \code{asinh} transformation.
#'
#'@param num_col An integer vector of index indicating the columns to be transform. 
#'Default is \code{1:ncol(M)} to transform all the data. 
#'
#'@return An object of class 'cytomeTree' providing a partitioning
#'of the set of n cells.
#'\itemize{
#'\item{\code{annotation}}{ A \code{data.frame} containing the annotation of each 
#'cell population underlying the tree pattern.}
#'\item{\code{labels}}{ The partitioning of the set of n cells.}
#'\item{\code{M}}{ The transformed matrix of mass cytometry.}
#'\item{\code{mark_tree}}{ A two level list containing markers used
#'for node splitting.}
#'\item{\code{transformation}}{ Transformation used}
#'\item{\code{num_col}}{ Indexes of columns transformed}
#' }
#'
#'@details First of all, data can be transformed using different transformations. 
#'The algorithm is based on the construction of a binary tree,
#'the nodes of which are subpopulations of cells. At each node,
#'observed cells and markers are modeled by both a family of normal
#'distributions and a family of bi-modal normal mixture distributions.
#'Splitting is done according to a normalized difference of AIC between
#'the two families.
#'@author Anthony Devaux, Boris Hejblum
#'
#'@export
#'
#'@examples
#'data(IMdata)
#'
#'# dimension of data
#'dim(IMdata)
#'
#'# given the size of the dataset, the code below can take several minutes to run
#'
#'\dontrun{
#'# Don't transform Time et Cell_length column
#'num_col <- 3:ncol(IMdata)
#'
#'# Build Cytoftree binary tree
#'tree <- CytofTree(M = IMdata, minleaf = 1, t = 0.1, transformation = "asinh", num_col = num_col)
#'
#'# Annotation
#'annot <- Annotation(tree, plot = FALSE, K2markers = colnames(IMdata))
#'
#'# Provide subpopulations
#'annot$combinations
#'}
#'

CytofTree <- function(M, minleaf = 1, t = .1, verbose = TRUE, force_first_markers = NULL, 
                      transformation = c("asinh", "biexp", "log10", "none"), 
                      num_col = 1:ncol(M)){
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
  
  if(any(!transformation%in%c("asinh", "biexp", "log10", "none"))){
    stop("The only allowed values for transformation are 'asinh', 'log10' or 'none'")
  }
  
  if(length(transformation)>1){
    transformation <- transformation[1]
    warning("Default 'asinh' transformation used")
  }
  
  if(transformation!="none"){
    if(any(!(num_col%in%1:ncol(M)))){
      stop(paste0("'num_col' should contain columns indices from 1 to ", ncol(M)))
    }
  }else{
    num_col <- NULL
  }
  
  if(transformation=="asinh"){
    M <- AsinhTransformation(M, num_col)
  }
  
  if(transformation=="biexp"){
    M <- BiexpTransformation(M, num_col)
  }
  
  if(transformation=="log10"){
    M <- Log10Transformation(M, num_col)
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
