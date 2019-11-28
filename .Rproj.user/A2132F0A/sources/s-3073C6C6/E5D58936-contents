#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::mat cpp_dmvnorm_invmat(arma::vec evals){
  int n = evals.n_elem;
  arma::mat output(n,n,fill::zeros);
  for (int i=0;i<n;i++){
    output(i,i) = 1/evals(i);
  }
  return(output);
}

// [[Rcpp::export]]
arma::vec cpp_dmvnorm(arma::mat &X, arma::vec mean, arma::mat Sigma){
  // parameters
  int N = X.n_rows;
  int K = X.n_cols;
  
  // eigenvalue decomposition
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, Sigma);
  
  double    detSig  = arma::prod(eigval);
  arma::mat diaginv = cpp_dmvnorm_invmat(eigval);
  arma::mat invSig  = eigvec*diaginv*eigvec.t();
  
  // compute with iteration
  double term1 = -0.5*std::log(detSig);
  double term2 = -(static_cast<double>(K)/2.0)*std::log(2.0*3.14159265358979323846);
  double term3 = term1+term2;
  arma::vec logout(N,fill::zeros);
  for (int n=0;n<N;n++){
    arma::vec xmu = X.row(n).t()-mean;
    logout(n) = arma::dot(invSig*xmu, xmu)*(-0.5) + term3;
  }
  
  // return the output
  return(logout);
}