#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector nmixpdfC(NumericVector x,
NumericVector mu,
NumericVector Sigma,
NumericVector Pi,
bool Log = true)
{
  int N = x.size();
  int K = mu.size();
  NumericVector y = NumericVector(N);
  for (int n = 0; n < N; n++){
    double temp = 0;
    for (int k = 0; k < K; k++)
    {
      double m = mu(k);
      double s = Sigma(k);
      double p = Pi(k);
      double cons = 1/sqrt(2.0 * M_PI * s);
      double comp = cons * exp(-.5 * ((x(n)-m)*(x(n)-m))/s);
      temp = temp + p * comp;
    }
    y(n) = temp;
  }
  if(Log)
  {
    y = log(y); 
  }
  return y;
}


