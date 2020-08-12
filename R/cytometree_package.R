#'@useDynLib cytometree, .registration = TRUE
#'@importFrom Rcpp sourceCpp evalCpp
#'
#'@details 
#'The main function in this package is \code{\link{CytomeTree}}.
#'\tabular{ll}{
#'Package: \tab cytometree\cr
#'Type: \tab Package\cr
#'Version: \tab 2.0.4\cr
#'Date: \tab 2020-08-12\cr
#'License:\tab \href{http://www.gnu.org/licenses/lgpl.txt}{LGPL-3}\cr
#'}
#'The algorithm is based on the construction of a binary tree, 
#'the nodes of which are subpopulations of cells. At each node, 
#'observed cells and markers are modeled by both a family of normal
#'distributions and a family of bi-modal normal mixture distributions. 
#'Splitting is done according to a normalized difference of AIC between 
#'the two families.
#'Given the unsupervised nature of the binary tree, some of the available
#'markers may not be used to find the different cell populations present in 
#'a given sample. To recover a complete annotation, we defined, as a post 
#'processing procedure, an annotation method which allows the user to 
#'distinguish two or three expression levels per marker.
#'
#'@references Commenges D, Alkhassim C, Gottardo R, Hejblum BP, Thiébaut R (2018). 
#'cytometree: a binary tree algorithm for automatic gating in cytometry analysis. 
#'Cytometry Part A, 93(11):1132-1140. <doi: 10.1002/cyto.a.23601>
#'
"_PACKAGE"
