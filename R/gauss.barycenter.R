#' Barycenter of Gaussian Distributions
#' 
#' @examples 
#' ## generate univariate Gaussians
#' mylist1d = list()
#' for (i in 1:50){
#'    my.ctd = stats::runif(1, min=-0.5, max=0.5)
#'    my.sd  = stats::runif(1, min=0.9, max=1.1)
#'    mylist1d[[i]] = wrapgauss1d(mean=my.ctd, sd=my.sd)
#' }
#' wass2fpt = gauss.barycenter(mylist1d, type="wass2fpt")
#' wass2rgd = gauss.barycenter(mylist1d, type="wass2rgd")
#' 
#' \dontrun{
#' ## test with 5-dimensional Gaussians
#' mylist5d = list()
#' mycovs   = list()
#' for (i in 1:50){
#'    my.ctd = stats::runif(5, min=-0.5, max=0.5)
#'    my.sig = stats::cov(matrix(T4Gauss::rnorm(100*5),ncol=5))
#'    mycovs[[i]] = my.sig
#'    mylist5d[[i]] = wrapgaussNd(mu=my.ctd, sigma=my.sig)
#' }
#' wass2fpt = gauss.barycenter(mylist5d, type="wass2fpt")
#' wass2rgd = gauss.barycenter(mylist5d, type="wass2rgd")
#' }
#' 
#' @export
gauss.barycenter <- function(x, y=NULL, type=c("wass2fpt","wass2rgd"), 
                             maxiter=200, eps=1e-10, nthreads=0){
  #######################################################
  # Preprocessing
  mymethod = match.arg(type)
  if (is.list(x)){
    if (!check_list_gauss(x)){
      stop("* gauss.barycenter : input 'x' should be a list of 'wrapgauss' objects having same dimension.")
    }
    dglist = x
  } else {
    xcond = base::inherits(x, "wrapgauss")
    ycond = base::inherits(y, "wrapgauss")
    if (!(xcond&&ycond)){
      stop(" gauss.barycenter : input 'x' and 'y' should be of class 'wrapgauss'.")
    }
    dglist = list()
    dglist[[1]] = x
    dglist[[2]] = y
    if (!check_list_gauss(dglist)){
      stop("* gauss.barycenter : 'x' and 'y' seem to be inconsistent.")
    }
  }
  par.iter = round(maxiter)
  par.eps  = as.double(eps)
  nCores   = round(nthreads)
  
  #######################################################
  # Set the weight
  lambdas = rep(1, length(dglist))/length(dglist)
  
  #######################################################
  # Compute and Return
  return(barygauss_selection(dglist, lambdas, mymethod, par.iter, par.eps, nCores))
}



# auxiliary functions -----------------------------------------------------
# (0) switching argument
# (1) wass2    : under Wasserstein Geometry
#     wass2fpt : fixed-point iteration
#     wass2rgd : Riemannian Optimization

# (0) switching argument --------------------------------------------------
#' @keywords internal
#' @noRd
barygauss_selection <- function(dglist, lbd=(rep(1,length(dglist))/length(dglist)), method, 
                                par.iter, par.eps, nCores=0){
  output = switch(method,
                  wass2fpt = barygauss_wass2fpt(dglist, lbd, par.iter, par.eps, nCores),
                  wass2rgd = barygauss_wass2rgd(dglist, lbd, par.iter, par.eps, nCores))
  return(output)
}

# (1) wass2 ---------------------------------------------------------------
#' @keywords internal
#' @noRd
barygauss_wass2fpt <- function(dglist, lambdas, par.iter, par.eps, nCores=0){
  # parameters
  N = length(dglist)
  p = dglist[[1]]$dimension
  
  # stack all means and covs : be careful about the dimension
  mean3 = array(0,c(N,p))
  covs3 = array(0,c(p,p,N))
  for (n in 1:N){
    mean3[n,]  = as.vector(dglist[[n]]$mu)
    covs3[,,n] = as.matrix(dglist[[n]]$sigma)
  }
  
  # compute 
  tmpout = wass2_barycenter(mean3, covs3, lambdas, par.iter, par.eps)
  if (p < 2){
    return(wrapgauss1d(mean=tmpout$mu_bary, sd=sqrt(tmpout$sig_bary)))
  } else {
    return(wrapgaussNd(mu=tmpout$mu_bary, sigma=tmpout$sig_bary))
  }
}
#' @keywords internal
#' @noRd
barygauss_wass2rgd <- function(dglist, lambdas, par.iter, par.eps, nCores=0){
  # parameters
  N = length(dglist)
  p = dglist[[1]]$dimension
  
  # stack all means and covs : be careful about the dimension
  mean3 = array(0,c(N,p))
  covs3 = array(0,c(p,p,N))
  for (n in 1:N){
    mean3[n,]  = as.vector(dglist[[n]]$mu)
    covs3[,,n] = as.matrix(dglist[[n]]$sigma)
  }
  
  # mean
  mout = rep(0,p)
  for (n in 1:N){
    mout = mout + lambdas[n]*as.vector(mean3[n,])
  }
  mout = as.vector(wass2covs_mu(mean3, lambdas))
  
  # covariance
  if (nCores > 0){
    covout = wass2covs_fmean_openmp(covs3, lambdas, par.iter, par.eps, nCores)$x  
  } else {
    covout = wass2covs_fmean(covs3, lambdas, par.iter, par.eps)$x
  }
  
  if (p < 2){
    return(wrapgauss1d(mean=mout, sd=sqrt(covout)))
  } else {
    return(wrapgaussNd(mu=mout, sigma=covout))
  }
}

# ## personal experiment with microbenchmark
# library(microbenchmark)
# mylist5d = list()
# mycovs   = list()
# for (i in 1:50){
#   my.ctd = stats::runif(5, min=-0.5, max=0.5)
#   my.sig = stats::cov(matrix(T4Gauss::rnorm(100*5),ncol=5))
#   mycovs[[i]] = my.sig
#   mylist5d[[i]] = wrapgaussNd(mu=my.ctd, sigma=my.sig)
# }
# microbenchmark(
#   fpt = gauss.barycenter(mylist5d, type="wass2fpt"),
#   rgd0 = gauss.barycenter(mylist5d, type="wass2rgd"),  
#   rgd3 = gauss.barycenter(mylist5d, type="wass2rgd", nthreads = 3),
#   rgd6 = gauss.barycenter(mylist5d, type="wass2rgd", nthreads = 6),
#   rgd9 = gauss.barycenter(mylist5d, type="wass2rgd", nthreads = 9)
# )