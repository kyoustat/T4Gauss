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
arma::vec wass2covs_normvec(arma::mat mold, arma::cube cov3, int nCores){
  // parameter
  int p = mold.n_rows;
  int N = cov3.n_slices;

  arma::mat tmp(p,p,fill::zeros);
  arma::vec output(N,fill::zeros);
  
  // iterate
  for (int n=0;n<N;n++){
    tmp = cov3.slice(n);
    output(n) = std::sqrt(static_cast<float>(arma::trace(mold + tmp - (2.0*arma::real(arma::sqrtmat(mold*tmp))))));
  }
  return(output);
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

// (2) frechet median of covariances under 2-Wasserstein geometry
// [[Rcpp::export]]
Rcpp::List wass2covs_fmedian(arma::cube covs, arma::vec lambdas, int maxiter, double eps){
  // get parameters
  int N = covs.n_slices;
  int p = covs.n_cols;
  int iter = 0;
  
  // initialize
  arma::mat mold = arma::mean(covs, 2);
  arma::mat mnew(p,p,fill::zeros);
  
  arma::cube tvecs; tvecs.copy_size(covs); tvecs.fill(0.0); 
  arma::mat  dtmp;  dtmp.copy_size(mold);  dtmp.fill(0.0); // on TpM
  arma::vec  normvec(N,fill::zeros);
  arma::uvec nonsingular;
  
  // let's iterate
  double sqnorm = 10000.00;
  while (sqnorm > eps){
    // 1-1. compute log-pulled vectors and norm
    tvecs.fill(0.0);
    for (int n=0;n<N;n++){
      tvecs.slice(n) = wass2covs_log(mold, covs.slice(n));
    }
    normvec = wass2covs_normvec(mold, covs, 0);
    // for (int n=0;n<N;n++){
    //   normvec(n) = (std::sqrt(static_cast<float>(arma::trace(tvecs.slice(n)*mold*tvecs.slice(n)))));
    // }
    
    // 1-2. find the one with non-singular distance
    nonsingular = arma::find(normvec>1e-10);
    int M = nonsingular.n_elem;
    if (M==0){
      break;
    }
    // 2-1. update numerator
    dtmp.fill(0.0);
    for (int j=0;j<M;j++){
      dtmp += tvecs.slice(nonsingular(j))*lambdas(nonsingular(j))/normvec(nonsingular(j));
    }
    // 2-2. update denominator
    double denom = 0.0;
    for (int j=0;j<M;j++){
      denom += lambdas(nonsingular(j))/normvec(nonsingular(j));
    }
    dtmp /= denom;
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
Rcpp::List wass2covs_fmedian_openmp(arma::cube covs, arma::vec lambdas, int maxiter, double eps, int nCores){
  // get parameters
  int N = covs.n_slices;
  int p = covs.n_cols;
  int iter = 0;
  
  // initialize
  arma::mat mold = arma::mean(covs, 2);
  arma::mat mnew(p,p,fill::zeros);
  
  arma::cube tvecs; tvecs.copy_size(covs); tvecs.fill(0.0); 
  arma::mat  dtmp;  dtmp.copy_size(mold);  dtmp.fill(0.0); // on TpM
  arma::vec  normvec(N,fill::zeros);
  arma::uvec nonsingular;
  
  // let's iterate
  double sqnorm = 10000.00;
  while (sqnorm > eps){
    // 1-1. compute log-pulled vectors and norm
    tvecs.fill(0.0);
    #pragma omp parallel for num_threads(nCores) shared(tvecs, mold, covs)
    for (int n=0;n<N;n++){
      tvecs.slice(n) = wass2covs_log(mold, covs.slice(n));
    }
    normvec = wass2covs_normvec(mold, covs, nCores);
    // for (int n=0;n<N;n++){
    //   normvec(n) = (std::sqrt(static_cast<float>(arma::trace(tvecs.slice(n)*mold*tvecs.slice(n)))));
    // }
    
    // 1-2. find the one with non-singular distance
    nonsingular = arma::find(normvec>1e-10);
    int M = nonsingular.n_elem;
    if (M==0){
      break;
    }
    // 2-1. update numerator
    dtmp.fill(0.0);
    for (int j=0;j<M;j++){
      dtmp += tvecs.slice(nonsingular(j))*lambdas(nonsingular(j))/normvec(nonsingular(j));
    }
    // 2-2. update denominator
    double denom = 0.0;
    for (int j=0;j<M;j++){
      denom += lambdas(nonsingular(j))/normvec(nonsingular(j));
    }
    dtmp /= denom;
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