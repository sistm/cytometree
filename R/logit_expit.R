#'logit and expit functions
#'
#'@keywords internal
#'

logit <- function(x){
  log(x/(1 - x))
}

expit <- function(x){
  1/(1 + exp(-x))
}
