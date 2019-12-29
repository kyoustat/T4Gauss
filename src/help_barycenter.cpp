/*
 * (1) wass2_barycenter
 */

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// (1) wass2_barycenter
// [[Rcpp::export]]
Rcpp::List wass2_barycenter(arma::mat mean3, arma::cube covs3, arma::vec lambdas, int maxiter, double eps){
  // parameters
  int N = mean3.n_rows;
  int p = mean3.n_cols;
  
  // compute 1. mean vector
  arma::vec mvec(p,fill::zeros);
  for (int n=0;n<N;n++){
    mvec += lambdas(n)*(mean3.row(n).t());
  }
  
  // compute 2. covariance matrix
  arma::mat msigold(p,p,fill::eye);
  arma::mat msignew(p,p,fill::zeros);
  arma::mat msohalf(p,p,fill::zeros); // Sig^(1/2)
  arma::mat msohinv(p,p,fill::zeros); // Sig^(-1/2)
  arma::mat msumtmp(p,p,fill::zeros); // intermediate term for addition
  
  double sthr = eps;
  double sinv = 0.0;
  for (int it=0;it<maxiter;it++){
    // precompute two interim matrices
    msohalf = arma::sqrtmat_sympd(msigold);
    msohinv = arma::inv_sympd(msohalf);
    
    // intermediate term for addition
    msumtmp.fill(0.0);
    for (int n=0;n<N;n++){
      msumtmp += lambdas(n)*arma::sqrtmat_sympd(msohalf*covs3.slice(n)*msohalf);
    }
    msignew = msohinv*msumtmp*msumtmp*msohinv;
    
    // updating rule
    sinv    = arma::norm(msignew-msigold,"fro");
    msigold = msignew;
    if (sinv < sthr){
      break;
    }
  }
  
  // return
  return Rcpp::List::create(Rcpp::Named("mu_bary")=mvec,
                            Rcpp::Named("sig_bary")=msigold);
}