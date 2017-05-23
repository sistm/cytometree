#'Binary tree algorithm for cytometry data analysis.
#'
#'The algorithm is based on the construction of a binary tree, 
#'the nodes of which are subpopulations of cells. At each node, 
#'observed cells and markers are modeled by both a family of normal
#'and a family of normal mixtures distributions. 
#'Splitting is done according to a normalized difference of AIC between 
#'the two families.
#'Given the unsupervised nature of the binary tree, some of the available
#'markers may not be used to find the different cell populations present in 
#'a given sample. To recover a complete annotation, we defined, as a post 
#'processing procedure, an annotation method which allows the user to 
#'distinguish two or three expression levels per marker.
#'
#'\tabular{ll}{
#'Package: \tab cytometree\cr
#'Type: \tab Package\cr
#'Version: \tab 1.0.0\cr
#'Date: \tab 2017-05-23\cr
#'License:\tab \href{http://www.gnu.org/licenses/lgpl.txt}{LGPL-3}\cr
#'}
#'The main function in this package is \code{\link{CytomeTree}}.
#'
#'@author Chariff Alkhassim, Boris P. Hejblum
#'--- Maintainer: Chariff Alkhassim
#'
#'@references D. Commenges, C Alkhassim, B.P. Hejblum, R. Thiébaut. A Binary tree algorithm for automatic gating in cytometry analysis. Sumitted, 2017.
#'
#'@docType package
#'@name cytometree-package
#'@aliases cytometree
#'
#'@useDynLib cytometree, .registration = TRUE
#'@importFrom Rcpp evalCpp
#'
NULL