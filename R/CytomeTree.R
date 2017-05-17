#' Binary tree algorithm for cytometry data analysis.
#'
#'@param M A matrix of size n x p containing cytometry measures
#'of n cells on p markers.  
#'
#'@param minleaf An integer indicating the minimum number of cell
#'per population. Default is \code{1}.
#'
#'@param t A real positive-or-null number used for comparison with
#'the normalized AIC computed at each node of the tree.
#'A higher values limits the height of the tree. 
#'
#'@return An object of class 'cytomeTree' providing a partitioning
#'of the set of n cells. 
#'\itemize{
#'\item{\code{combinations}}{ A matrix of size n x p containing 
#'the annotation of each cell given by the tree.}
#'\item{\code{labels}}{ The partitioning of the set of n cells.}
#'\item{\code{M}}{ The input matrix.}
#'\item{\code{mark_tree}}{ A two level list containing markers used 
#'for node splitting.}
#' }
#'
#'@details The algorithm is based on the construction of a binary tree, 
#'the nodes of which are subpopulations of cells. At each node, 
#'observed cells and markers are modeled by both a familly of normal
#'and a familly of normal mixtures distributions. 
#'Spliting is done according to a normalized difference of AIC between 
#'the two families.
#'@author Chariff Alkhassim
#'
#'@export 
# 
CytomeTree <- function(M, minleaf = 1, t = .1)
{
  if((class(M) != "matrix") & (class(M) != "data.frame"))
  {
    stop("M should be of class matrix or data.frame")
  }
  n <- nrow(M)
  if(minleaf >= n)
  {
    stop("minleaf is superior to n")
  }
  p <- ncol(M)
  if(p > n){
    stop("p is superior to n")
  }
  if(any(is.na(M)))
  {
    stop("M contains NAs")
  }
  BT <- BinaryTree(M, minleaf, t)
  Tree <- list("M" = M, "labels" = BT$labels,
               "pl_list"= BT$pl_list,
               "mark_tree" = BT$mark_tree,
               "combinations" = BT$combinations)
  class(Tree) <- "CytomeTree"
  return(Tree)
}