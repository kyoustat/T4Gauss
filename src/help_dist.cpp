/*
 * (1) wass2_dist  : pairwise distance of a set
 * (2) wass2_dist2 : pairwise distance of two sets 
 */

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// (1) wass2_dist
// [[Rcpp::export]]
arma::mat wass2_dist(arma::mat mean3, arma::cube covs3){
  // parameters
  int N = mean3.n_rows;
  int p = mean3.n_cols;
  
  // precompute sqrtm
  arma::cube sqrt3(p,p,N,fill::zeros);
  for (int n=0;n<N;n++){
    sqrt3.slice(n) = arma::sqrtmat_sympd(covs3.slice(n));
  }
  
  // main computation
  arma::rowvec m0(p,fill::zeros);
  arma::rowvec m1(p,fill::zeros);
  arma::mat sig0(p,p,fill::zeros);
  arma::mat sig1(p,p,fill::zeros);
  arma::mat sig0half(p,p,fill::zeros);
  
  arma::mat output(N,N,fill::zeros);
  for (int n=0;n<(N-1);n++){
    m0   = mean3.row(n);
    sig0 = covs3.slice(n);
    sig0half = sqrt3.slice(n);
    for (int m=(n+1);m<N;m++){
      m1   = mean3.row(m);
      sig1 = covs3.slice(m);
      output(n,m) = std::pow(arma::norm(m0-m1,2),2) + arma::trace(sig0 + sig1 - 2.0*arma::sqrtmat_sympd(sig0half*sig1*sig0half));
      output(m,n) = output(n,m);
    }
  }
  
  // return
  return(output);
}



// (2) wass2_dist2
// [[Rcpp::export]]
arma::mat wass2_dist2(arma::mat mean1, arma::mat mean2, arma::cube covs1, arma::cube covs2){
  // parameters
  int M = mean1.n_rows;
  int N = mean2.n_rows;
  int p = mean1.n_cols;
  
  // precompute sqrtm of 1st set of covariances
  arma::cube sqrt3(p,p,M,fill::zeros);
  for (int m=0;m<M;m++){
    sqrt3.slice(m) = arma::sqrtmat_sympd(covs1.slice(m));
  }
  
  // main computation
  arma::rowvec m0(p,fill::zeros);
  arma::rowvec m1(p,fill::zeros);
  arma::mat sig0(p,p,fill::zeros);
  arma::mat sig1(p,p,fill::zeros);
  arma::mat sig0half(p,p,fill::zeros);
  
  arma::mat output(M,N,fill::zeros);
  for (int m=0;m<M;m++){
    m0   = mean1.row(m);
    sig0 = covs1.slice(m);
    sig0half = sqrt3.slice(m);
    for (int n=0;n<N;n++){
      m1   = mean2.row(n);
      sig1 = covs2.slice(n);
      output(m,n) = std::pow(arma::norm(m0-m1,2),2) + arma::trace(sig0 + sig1 - 2.0*arma::sqrtmat_sympd(sig0half*sig1*sig0half));
    }
  }
  
  // return
  return(output);
}
