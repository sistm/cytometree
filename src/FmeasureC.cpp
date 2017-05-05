#include <RcppArmadillo.h>
// #define ARMA_64BIT_WORD // to enable matrix with more than 4 billions elements
// requires a 64bits machine
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' C++ implementation of the F-measure computation
//' 
//'@param pred vector of a predicted partition
//'@param ref vector of a reference partition
//'
//'@export
//'
//'@examples
//'c <- list(c(1,1,2,3,2,3), c(1,1,1,2,3,3),c(2,2,1,1,1,1))
//'similarityMatC(sapply(c, "["))
//'
//'c2 <- list()
//'for(i in 1:100){
//'     c2 <- c(c2, list(rmultinom(n=1, size=3000, prob=rexp(n=3000))))
//'}
//'library(microbenchmark)
//'f <- function(){c3 <-sapply(c2, "[")
//'             similarityMatC(c3)}
//'microbenchmark(f(), time=1L)
//'
// [[Rcpp::export]]
double FmeasureC(NumericVector pred, NumericVector ref){
  
  vec K = unique(pred);
  K = sort(K);
  vec C = unique(ref);
  C = sort(C);
  int m = K.size();
  int n = C.size();
  
  mat M = mat(n, m);
  mat Pr = mat(n, m);
  mat Re = mat(n, m);
  mat Fmat = mat(n, m);
  
  vec C_card = vec(n);
  vec K_card = vec(m);
  
  for(int i=0; i<n; i++){
    C_card(i) = sum(ref == C(i));
    for(int j=0; j<m; j++){
      K_card(j) = sum(pred == K(j));
      M(i,j) = sum((ref==C(i)) & (pred==K(j)));
      Pr(i,j) = M(i,j)/K_card(j);
      Re(i,j) = M(i,j)/C_card(i);
      if((Pr(i,j) + Re(i,j)) == 0.0){
        Fmat(i,j) = 0;
      }else{
        Fmat(i,j) = 2.0*Pr(i,j)*Re(i,j)/(Pr(i,j) + Re(i,j));
      }
    }
  } 
  
  double C_card_sum = sum(C_card);
  vec Ffinal = vec(n);
  vec Fsum = vec(n);
  
  for(int i=0; i<n; i++){
    Ffinal(i) = max(Fmat.row(i));
    Fsum(i) = Ffinal(i)*C_card(i)/C_card_sum;
  }
  double Ftotal = sum(Fsum);
  
  return Ftotal;
}


//' C++ implementation of the F-measure computation without the ref classe 0
//' 
//' Aghaeepour in FlowCAP 1 ignore the reference class labeled "0"
//' 
//' 
//'@param pred vector of a predicted partition
//'@param ref vector of a reference partition
//'
//'@export
//'
//'@examples
//'c <- list(c(1,1,2,3,2,3), c(1,1,1,2,3,3),c(2,2,1,1,1,1))
//'similarityMatC(sapply(c, "["))
//'
//'c2 <- list()
//'for(i in 1:100){
//'     c2 <- c(c2, list(rmultinom(n=1, size=3000, prob=rexp(n=3000))))
//'}
//'library(microbenchmark)
//'f <- function(){c3 <-sapply(c2, "[")
//'             similarityMatC(c3)}
//'microbenchmark(f(), time=1L)
//'
// [[Rcpp::export]]
double FmeasureC_no0(NumericVector pred, NumericVector ref){
  
  vec K = unique(pred);
  K = sort(K);
  vec C = unique(ref);
  C = sort(C);
  int p = K.size();
  int n = C.size();
  vec C_no0 = C.subvec(1, n-1);
  int n_no0 = n-1;
  
  mat M = mat(n_no0, p);
  mat Pr = mat(n_no0, p);
  mat Re = mat(n_no0, p);
  mat Fmat = mat(n_no0, p);
  
  vec C_card = vec(n_no0);
  vec K_card = vec(p);
  
  for(int i=0; i<n_no0; i++){
    C_card(i) = sum(ref == C_no0(i));
    for(int j=0; j<p; j++){
      K_card(j) = sum((pred == K(j)) & (ref !=C(0)));
      M(i,j) = sum((ref==C_no0(i)) & (pred==K(j)));
      Pr(i,j) = M(i,j)/K_card(j);
      Re(i,j) = M(i,j)/C_card(i);
      if((Pr(i,j) + Re(i,j)) == 0.0){
        Fmat(i,j) = 0;
      }else{
        Fmat(i,j) = 2.0*Pr(i,j)*Re(i,j)/(Pr(i,j) + Re(i,j));
      }
    }
  } 
  
  double C_card_sum = sum(C_card);
  vec Ffinal = vec(n_no0);
  vec Fsum = vec(n_no0);
  
  for(int i=0; i<n_no0; i++){
    Ffinal(i) = max(Fmat.row(i));
    Fsum(i) = Ffinal(i)*C_card(i)/C_card_sum;
  }
  double Ftotal = sum(Fsum);
  
  return Ftotal;
}


//' C++ implementation of cost computation with Fmeasure as loss function
//' 
//'
//'@param c matrix where each column is one MCMC partition
//'
//'@export
//'
//'@examples
//'c <- list(c(1,1,2,3,2,3), c(1,1,1,2,3,3),c(2,2,1,1,1,1))
//'Fmeasure_costC(sapply(c, "["))
//'
//'c2 <- list()
//'for(i in 1:500){
//'     c2 <- c(c2, list(rmultinom(n=1, size=10000, prob=rexp(n=10000))))
//'}
//'library(microbenchmark)
//'f <- function(){c3 <-sapply(c2, "[")
//'             Fmeasure_costC(c3)}
//'fa <- function(){c3 <-sapply(c2, "[")
//'             Fmeasure_costC_arma(c3)}
//'microbenchmark(f(), fa(), times=10L)
//'
// [[Rcpp::export]]
List Fmeasure_costC(NumericMatrix c){
  
  //int N = cc.n_cols;
  //int N = cc.n_rows;
  int N = c.ncol();
  //int n = c.nrow();
  
  
  NumericVector cost = NumericVector(N);
  mat Fmeas = mat(N, N, fill::eye);
  
  for(int i=0; i<N-1; i++){
    for(int j=i+1; j<N; j++){
      NumericVector pred_i = c(_,i);
      NumericVector ref_j = c(_,j);
      Fmeas(i,j) = FmeasureC(pred_i, ref_j);
      Fmeas(j,i) = FmeasureC(ref_j, pred_i); // WATCHOUT: FmeasureC is not symetric
    }
  }
  for(int k=0; k<N; k++){
    cost(k) = 1-(sum(Fmeas.col(k))-1)/N;
  }
  return Rcpp::List::create(Rcpp::Named("Fmeas") = Fmeas,
  Rcpp::Named("cost")=cost);
}

//' C++ implementation of cost computation with Fmeasure as loss function
//' using the Armadillo library
//' 
//'
//'@param c matrix where each column is one partition
//'
//'@export
// [[Rcpp::export]]
List Fmeasure_costC_arma(arma::mat c){
  
  int N = c.n_cols;
  //int n = c.n_rows;
  
  
  NumericVector cost = NumericVector(N);
  mat Fmeas = mat(N, N, fill::eye);
  
  for(int i=0; i<N-1; i++){
    for(int j=i+1; j<N; j++){
      vec pred_i_temp = c.col(i);
      vec ref_j_temp = c.col(j);
      NumericVector pred_i = as<NumericVector>(wrap(pred_i_temp));
      NumericVector ref_j = as<NumericVector>(wrap(ref_j_temp));
      Fmeas(i,j) = FmeasureC(pred_i, ref_j);
      Fmeas(j,i) = Fmeas(i,j);
    }
  }
  for(int k=0; k<N; k++){
    cost(k) = 1-(sum(Fmeas.col(k))-1)/N;
  }
  return Rcpp::List::create(Rcpp::Named("Fmeas") = Fmeas,
  Rcpp::Named("cost")=cost);
}


