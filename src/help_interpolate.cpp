/*
 * (1) wass2_interp
 */

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// (1) wass2_interp
// [[Rcpp::export]]
Rcpp::List wass2_interp(arma::vec m0, arma::vec m1, arma::mat sig0, arma::mat sig1, double t){
  // parameters
  int p = m0.n_elem;
  
  // 1. mean
  arma::vec mu_t = m0*(1.0-t) + t*m1;
  
  // 2. covariance
  arma::mat sig0half = arma::sqrtmat_sympd(sig0);
  arma::mat sig0hinv = arma::inv_sympd(sig0half);
  
  arma::mat sigmid = ((1.0-t)*sig0 + t*arma::sqrtmat_sympd(sig0half*sig1*sig0half));
  arma::mat sig_t  = sig0hinv*sigmid*sigmid*sig0hinv;
  
  // return output
  
  return Rcpp::List::create(Rcpp::Named("mu_t")=mu_t,
                            Rcpp::Named("sig_t")=sig_t);
}
