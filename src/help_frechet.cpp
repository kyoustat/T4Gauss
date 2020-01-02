/*
 * (1) frechet mean   of covariances under 2-Wasserstein geometry
 * (2) frechet median of covariances under 2-Wasserstein geometry
 */

#include "RcppArmadillo.h"
#include <omp.h>

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// auxiliary functions
arma::mat wass2covs_lyapunov(arma::mat bot, arma::mat upp){
  arma::mat output = arma::syl(bot,bot,-upp);
  return(output);
}
arma::mat wass2covs_log(arma::mat C, arma::mat B){ // C a center, B another point
  return(arma::real(arma::sqrtmat(B*C)) + arma::real(arma::sqrtmat(C*B)) - 2.0*C);
}
arma::mat wass2covs_exp(arma::mat C, arma::mat V){
  int p = C.n_rows;
  arma::mat diagp(p,p,fill::eye);
  arma::mat LCV = wass2covs_lyapunov(C,V);
  return((LCV+diagp)*C*(LCV+diagp));
}

// (1) frechet mean of covariances under 2-Wasserstein geometry
// [[Rcpp::export]]
arma::rowvec wass2covs_mu(arma::mat means, arma::vec lambdas){
  int N = means.n_rows;
  int p = means.n_cols;
  
  arma::rowvec output(p,fill::zeros);
  for (int n=0;n<N;n++){
    output += means.row(n)*lambdas(n);
  }
  return(output);
}
// [[Rcpp::export]]
Rcpp::List wass2covs_fmean(arma::cube covs, arma::vec lambdas, int maxiter, double eps){
  // get parameters
  int N = covs.n_slices;
  int p = covs.n_cols;
  int iter = 0;
  
  // initialize
  arma::mat mold = arma::mean(covs, 2);
  arma::mat mnew(p,p,fill::zeros);
  
  arma::cube tvecs; tvecs.copy_size(covs); tvecs.fill(0.0); 
  arma::mat  dtmp;  dtmp.copy_size(mold);  dtmp.fill(0.0); // on TpM
  
  // let's iterate
  double sqnorm = 10000.00;
  while (sqnorm > eps){
    // 1. compute log-pulled vectors
    tvecs.fill(0.0);
    for (int n=0;n<N;n++){
      tvecs.slice(n) = wass2covs_log(mold, covs.slice(n));
    }
    // 2. compute updating scheme with weighted version
    dtmp.fill(0.0);
    for (int n=0;n<N;n++){
      dtmp += tvecs.slice(n)*lambdas(n); 
    }
    sqnorm = arma::norm(dtmp, "fro");
    if (sqnorm > eps){
      // 3. update using exponential map
      mnew = wass2covs_exp(mold, dtmp);
      // 4. iteration : update sqnorm
      sqnorm = arma::norm(dtmp, "fro");
      // 6. update others
      mold = mnew;
    }
    // 5. iteration : iter
    iter += 1;
    if (iter >= maxiter){
      break;
    }
  }
  
  // return
  return(Rcpp::List::create(Rcpp::Named("x")=mold,
                            Rcpp::Named("iteration")=iter));
}
// [[Rcpp::export]]
Rcpp::List wass2covs_fmean_openmp(arma::cube covs, arma::vec lambdas, int maxiter, double eps, int nCores){
  // get parameters
  int N = covs.n_slices;
  int p = covs.n_cols;
  int iter = 0;
  
  // initialize
  arma::mat mold = arma::mean(covs, 2);
  arma::mat mnew(p,p,fill::zeros);
  
  arma::cube tvecs; tvecs.copy_size(covs); tvecs.fill(0.0); 
  arma::mat  dtmp;  dtmp.copy_size(mold);  dtmp.fill(0.0); // on TpM
  
  // let's iterate
  double sqnorm = 10000.00;
  while (sqnorm > eps){
    // 1. compute log-pulled vectors
    tvecs.fill(0.0);
    #pragma omp parallel for num_threads(nCores) shared(tvecs, mold, covs)
    for (int n=0;n<N;n++){
      tvecs.slice(n) = wass2covs_log(mold, covs.slice(n));
    }
    // 2. compute updating scheme with weighted version
    dtmp.fill(0.0);
    for (int n=0;n<N;n++){
      dtmp += tvecs.slice(n)*lambdas(n); 
    }
    sqnorm = arma::norm(dtmp, "fro");
    if (sqnorm > eps){
      // 3. update using exponential map
      mnew = wass2covs_exp(mold, dtmp);
      // 4. iteration : update sqnorm
      sqnorm = arma::norm(dtmp, "fro");
      // 6. update others
      mold = mnew;
    }
    // 5. iteration : iter
    iter += 1;
    if (iter >= maxiter){
      break;
    }
  }
  
  // return
  return(Rcpp::List::create(Rcpp::Named("x")=mold,
                            Rcpp::Named("iteration")=iter));
}