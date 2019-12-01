#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List arma_gmm_full(arma::mat &X, int k, int maxiter){ // data are stacked as columns
  arma::gmm_full model;
  bool status = model.learn(X, k, maha_dist, random_subset, 10, maxiter, 1e-11, false);
  if (status==false){
    Rcpp::stop("* Fitting GMM with full covairance failed.");
  } else {
    // successful, let's take out the elements and return
    return Rcpp::List::create(Rcpp::Named("means")=model.means,
                              Rcpp::Named("covs")=model.fcovs,
                              Rcpp::Named("weight")=model.hefts,
                              Rcpp::Named("loglkd")=model.sum_log_p(X));
  }
}

// [[Rcpp::export]]
Rcpp::List arma_gmm_diag(arma::mat &X, int k, int maxiter){ // data are stacked as columns
  arma::gmm_diag model;
  bool status = model.learn(X, k, maha_dist, random_subset, 10, maxiter, 1e-11, false);
  if (status==false){
    Rcpp::stop("* Fitting GMM with diagonal covairance failed.");
  } else {
    // successful, let's take out the elements and return
    return Rcpp::List::create(Rcpp::Named("means")=model.means,
                              Rcpp::Named("covs")=model.dcovs,
                              Rcpp::Named("weight")=model.hefts,
                              Rcpp::Named("loglkd")=model.sum_log_p(X));
  }
}