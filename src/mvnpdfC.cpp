#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector mvnpdfC(NumericVector x,
                      double mean,
                      double var,
                      bool Log=true){
  int n  = x.size() ;
  NumericVector y = NumericVector(n);
  double constant = 1/sqrt(2.0 * M_PI *var);
  double constant_log = log(2.0 * M_PI *var);
  if(Log){
    for (int i = 0; i < n; i++){
      y(i) =  -0.5*(constant_log + ((x(i)-mean)* (x(i)-mean))/var);
    }
  }
  else{
    for (int i = 0; i < n; i++){
      y(i) =  constant * exp(-.5 * ((x(i)-mean)*(x(i)-mean))/var);
    }
  }
  return y;
}