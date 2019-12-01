/*
 * (1) sqrtm_covs : compute sqrtm for stacked 3d covariance matrices
 * 
 */
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// (1) sqrtm_covs
// [[Rcpp::export]]
arma::cube sqrtm_covs(arma::cube &covs3){
  // paramters & get ready
  int p = covs3.n_rows;
  int N = covs3.n_slices;
  arma::cube output(p,p,N,fill::zeros);
  
  // iteration
  arma::mat cslice(p,p,fill::zeros);
  for (int n=0;n<N;n++){
    cslice = covs3.slice(n);
    output.slice(n) = arma::sqrtmat_sympd(cslice);
  }
  
  // return
  return(output);
}