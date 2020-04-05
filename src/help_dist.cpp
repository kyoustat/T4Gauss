/*
 * (1) wass2_dist   : 2-Wasserstein Distances
 *     wass2_dist2 
 * (2) kl_dist      : Kullbeck-Leibler Divergence
 *     kl_dist2 
 * (3) skl_dist2    : Symmetrized KL (KL(p|q) + KL(q|p))
 * (4) cs_dist      : Cauchy-Schwarz Divergence
 *     cs_dist2
 * (5) bh_dist      : Bhattacharya
 *     bh_dist2 
 */

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// (1) Wasserstein Distances 
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
// [[Rcpp::export]]
arma::mat wass2_dist2(arma::mat mean1, arma::mat mean2, arma::cube covs1, arma::cube covs2){
  // parameters
  int M = mean1.n_rows;
  int N = mean2.n_rows;
  int p = mean1.n_cols;
  
  // precompute sqrtm of 1st set of covariances
  arma::cube sqrt3(p,p,N,fill::zeros);
  for (int n=0;n<N;n++){
    sqrt3.slice(n) = arma::sqrtmat_sympd(covs2.slice(n));
  }
  
  // main computation
  arma::rowvec m1(p,fill::zeros);
  arma::rowvec m2(p,fill::zeros);
  arma::mat sig1(p,p,fill::zeros);
  arma::mat sig2(p,p,fill::zeros);
  arma::mat sig2half(p,p,fill::zeros);
  
  arma::mat output(M,N,fill::zeros);
  for (int m=0;m<M;m++){
    m1   = mean1.row(m);
    sig1 = covs1.slice(m);
    for (int n=0;n<N;n++){
      m2   = mean2.row(n);
      sig2 = covs2.slice(n);
      output(m,n) = std::pow(arma::norm(m1-m2,2),2) + arma::trace(sig1 + sig2 - 2.0*arma::sqrtmat_sympd(sig2half*sig1*sig2half));
    }
  }
  
  // return
  return(output);
}
// (2) Kullbeck-Leibler Divergence
// [[Rcpp::export]]
arma::mat kl_dist(arma::mat mean3, arma::cube covs3){
  // parameters
  int N = mean3.n_rows;
  int p = mean3.n_cols;
  
  // precompute : inverse
  arma::cube pre_invs(p,p,N,fill::zeros);
  for (int n=0;n<N;n++){
    pre_invs.slice(n) = arma::inv_sympd(covs3.slice(n));
  }
  // precompute : determinant
  arma::vec pre_dets(N,fill::zeros);
  for (int n=0;n<N;n++){
    pre_dets(n) = arma::det(covs3.slice(n));
  }
  
  // main computation
  arma::colvec m1(p,fill::zeros);
  arma::colvec m2(p,fill::zeros);
  arma::colvec mdiff(p,fill::zeros);
  
  arma::mat S1(p,p,fill::zeros);
  arma::mat S2inv(p,p,fill::zeros);
  arma::mat Ip(p,p,fill::eye);
  
  double term1=0.0;
  double term2=0.0;
  arma::mat output(N,N,fill::zeros);
  for (int n=0;n<N;n++){
    m1 = mean3.row(n).t();
    S1 = covs3.slice(n);
    for (int m=0;m<N;m++){
      if (n!=m){
        m2    = mean3.row(m).t();
        S2inv = pre_invs.slice(m);
        
        mdiff = m1-m2;
        term1 = arma::dot(mdiff, S2inv*mdiff)/2.0;
        term2 = (arma::trace(S2inv*S1 - Ip) + std::log(pre_dets(m)) - std::log(pre_dets(n)))/2.0;
        
        output(n,m) = term1+term2;
      }
    }
  }
  
  // return
  return(output);
}
// [[Rcpp::export]]
arma::mat kl_dist2(arma::mat mean1, arma::cube covs1, arma::mat mean2, arma::cube covs2){
  // parameters
  int M = mean1.n_rows;
  int N = mean2.n_rows;
  int p = mean1.n_cols;
  
  // precompute : inverse
  arma::cube pre_invs2(p,p,N,fill::zeros);
  for (int n=0;n<N;n++){
    pre_invs2.slice(n) = arma::inv_sympd(covs2.slice(n));
  }
  // precompute : determinant
  arma::vec pre_dets1(M,fill::zeros);
  arma::vec pre_dets2(N,fill::zeros);
  for (int m=0;m<M;m++){
    pre_dets1(m) = arma::det(covs1.slice(m));
  }
  for (int n=0;n<N;n++){
    pre_dets2(n) = arma::det(covs2.slice(n));
  }
  
  // main computation
  arma::colvec m1(p,fill::zeros);
  arma::colvec m2(p,fill::zeros);
  arma::colvec mdiff(p,fill::zeros);
  
  arma::mat S1(p,p,fill::zeros);
  arma::mat S2inv(p,p,fill::zeros);
  arma::mat Ip(p,p,fill::eye);
  
  double term1=0.0;
  double term2=0.0;
  arma::mat output(M,N,fill::zeros);
  for (int m=0;m<M;m++){
    m1 = mean1.row(m).t();
    S1 = covs1.slice(m);
    for (int n=0;n<N;n++){
      m2    = mean2.row(n).t();
      S2inv = pre_invs2.slice(n);
      
      mdiff = m1-m2;
      term1 = arma::dot(mdiff, S2inv*mdiff)/2.0;
      term2 = (arma::trace(S2inv*S1 - Ip) + std::log(pre_dets2(n)) - std::log(pre_dets1(m)))/2.0;
      
      output(m,n) = term1+term2;
    }
  }
  // return
  return(output);
}
// (3) Symmetrized KL
// [[Rcpp::export]]
arma::mat skl_dist2(arma::mat mean1, arma::cube covs1, arma::mat mean2, arma::cube covs2){
  // parameters
  int M = mean1.n_rows;
  int N = mean2.n_rows;
  int p = mean1.n_cols;
  
  // precompute : inverse
  arma::cube pre_invs1(p,p,M,fill::zeros);
  arma::cube pre_invs2(p,p,N,fill::zeros);
  for (int m=0;m<M;m++){
    pre_invs1.slice(m) = arma::inv_sympd(covs1.slice(m));
  }
  for (int n=0;n<N;n++){
    pre_invs2.slice(n) = arma::inv_sympd(covs2.slice(n));
  }
  // precompute : determinant
  arma::vec pre_dets1(M,fill::zeros);
  arma::vec pre_dets2(N,fill::zeros);
  for (int m=0;m<M;m++){
    pre_dets1(m) = arma::det(covs1.slice(m));
  }
  for (int n=0;n<N;n++){
    pre_dets2(n) = arma::det(covs2.slice(n));
  }
  
  // main computation
  arma::colvec m1(p,fill::zeros);
  arma::colvec m2(p,fill::zeros);
  arma::colvec mdiff(p,fill::zeros);
  
  arma::mat S1(p,p,fill::zeros);
  arma::mat S2(p,p,fill::zeros);
  arma::mat S1inv(p,p,fill::zeros);
  arma::mat S2inv(p,p,fill::zeros);
  arma::mat Ip(p,p,fill::eye);
  
  double term1=0.0;
  double term2=0.0;
  double term3=0.0;
  double term4=0.0;
  arma::mat output(M,N,fill::zeros);
  for (int m=0;m<M;m++){
    m1    = mean1.row(m).t();
    S1    = covs1.slice(m);
    S1inv = pre_invs1.slice(m);
    for (int n=0;n<N;n++){
      m2    = mean2.row(n).t();
      S2    = covs2.slice(n);
      S2inv = pre_invs2.slice(n);
      
      mdiff = m1-m2;
      term1 = arma::dot(mdiff, S2inv*mdiff)/2.0;
      term2 = (arma::trace(S2inv*S1 - Ip) + std::log(pre_dets2(n)) - std::log(pre_dets1(m)))/2.0;
      term3 = arma::dot(mdiff, S1inv*mdiff)/2.0;
      term4 = (arma::trace(S1inv*S2 - Ip) + std::log(pre_dets1(m)) - std::log(pre_dets2(n)))/2.0;
      
      output(m,n) = term1+term2+term3+term4;
    }
  }
  // return
  return(output);
}
// (4) Cauchy-Schwarz Divergence
double cs_common(arma::vec m1, arma::vec m2, arma::mat S1, arma::mat S2){
  int p = m1.n_elem;
  arma::vec mdiff=m1-m2;
  
  double term1 = 0.0;
  double term2 = 0.0;
  double term3 = 0.0;
  
  term1 = std::log(arma::det(4.0*S1*S2))/(-2.0);
  term2 = std::log(arma::det(S1+S2));
  term3 = arma::dot(mdiff, arma::solve(S1+S2, mdiff));
  
  return(term1+term2+term3);
}
// [[Rcpp::export]]
arma::mat cs_dist(arma::mat mean3, arma::cube covs3){
  // parameters
  int N = mean3.n_rows;
  int p = mean3.n_cols;
  
  // main computation
  arma::colvec m1(p,fill::zeros);
  arma::colvec m2(p,fill::zeros);
  arma::mat S1(p,p,fill::zeros);
  arma::mat S2(p,p,fill::zeros);
  
  arma::mat output(N,N,fill::zeros);
  for (int m=0;m<(N-1);m++){
    m1  = mean3.row(m).t();
    S1  = covs3.slice(m);
    for (int n=(m+1);n<N;n++){
      m2 = mean3.row(n).t();
      S2 = covs3.slice(n);
      
      output(n,m) = cs_common(m1,m2,S1,S2);
      output(m,n) = output(n,m);
    }
  }
  
  // return
  return(output);
}
// [[Rcpp::export]]
arma::mat cs_dist2(arma::mat mean1, arma::cube covs1, arma::mat mean2, arma::cube covs2){
  // parameters
  int M = mean1.n_rows;
  int N = mean2.n_rows;
  int p = mean1.n_cols;
  
  // main computation
  arma::colvec m1(p,fill::zeros);
  arma::colvec m2(p,fill::zeros);
  arma::mat S1(p,p,fill::zeros);
  arma::mat S2(p,p,fill::zeros);
  
  arma::mat output(M,N,fill::zeros);
  for (int m=0;m<M;m++){
    m1  = mean1.row(m).t();
    S1  = covs1.slice(m);
    for (int n=0;n<N;n++){
      m2 = mean2.row(n).t();
      S2 = covs2.slice(n);
      
      output(m,n) = cs_common(m1,m2,S1,S2);
    }
  }
  
  // return
  return(output);
}
// (5) BH for Bhattacharya
double bh_common(arma::vec m1, arma::vec m2, arma::mat S1, arma::mat S2){
  int p = m1.n_elem;
  arma::vec mdiff=m1-m2;
  arma::mat Shalf=(S1+S2)/2.0;
  
  double term1 = 0.0;
  double term2top = 0.0;
  double term2bot = 0.0;
  
  term1 = arma::dot(mdiff, arma::solve(Shalf,mdiff))/8.0;
  term2top = arma::det(Shalf);
  term2bot = std::sqrt(arma::det(S1)*arma::det(S2));
  
  double output = 0.0;
  output = term1 + 0.5*(std::log(term2top) - std::log(term2bot));
  return(output);
}
// [[Rcpp::export]]
arma::mat bh_dist(arma::mat mean3, arma::cube covs3){
  // parameters
  int N = mean3.n_rows;
  int p = mean3.n_cols;
  
  // main computation
  arma::colvec m1(p,fill::zeros);
  arma::colvec m2(p,fill::zeros);
  arma::mat S1(p,p,fill::zeros);
  arma::mat S2(p,p,fill::zeros);
  
  arma::mat output(N,N,fill::zeros);
  for (int m=0;m<(N-1);m++){
    m1  = mean3.row(m).t();
    S1  = covs3.slice(m);
    for (int n=(m+1);n<N;n++){
      m2 = mean3.row(n).t();
      S2 = covs3.slice(n);
      
      output(n,m) = bh_common(m1,m2,S1,S2);
      output(m,n) = output(n,m);
    }
  }
  
  // return
  return(output);
}
// [[Rcpp::export]]
arma::mat bh_dist2(arma::mat mean1, arma::cube covs1, arma::mat mean2, arma::cube covs2){
  // parameters
  int M = mean1.n_rows;
  int N = mean2.n_rows;
  int p = mean1.n_cols;
  
  // main computation
  arma::colvec m1(p,fill::zeros);
  arma::colvec m2(p,fill::zeros);
  arma::mat S1(p,p,fill::zeros);
  arma::mat S2(p,p,fill::zeros);
  
  arma::mat output(M,N,fill::zeros);
  for (int m=0;m<M;m++){
    m1  = mean1.row(m).t();
    S1  = covs1.slice(m);
    for (int n=0;n<N;n++){
      m2 = mean2.row(n).t();
      S2 = covs2.slice(n);
      
      output(m,n) = bh_common(m1,m2,S1,S2);
    }
  }
  
  // return
  return(output);
}