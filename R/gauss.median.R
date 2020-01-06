#' Geometric Median of Gaussian Distributions
#' 
#' 
#' @export
gauss.median <- function(x, y=NULL, type=c("wass2"), maxiter=200, eps=1e-10, nthreads=1){
  #######################################################
  # Preprocessing
  mymethod = match.arg(type)
  if (is.list(x)){
    if (!check_list_gauss(x)){
      stop("* gauss.median : input 'x' should be a list of 'wrapgauss' objects having same dimension.")
    }
    dglist = x
  } else {
    xcond = base::inherits(x, "wrapgauss")
    ycond = base::inherits(y, "wrapgauss")
    if (!(xcond&&ycond)){
      stop(" gauss.median : input 'x' and 'y' should be of class 'wrapgauss'.")
    }
    dglist = list()
    dglist[[1]] = x
    dglist[[2]] = y
    if (!check_list_gauss(dglist)){
      stop("* gauss.median : 'x' and 'y' seem to be inconsistent.")
    }
  }
  par.iter = round(maxiter)
  par.eps  = as.double(eps)
  nCores   = round(nthreads)
  
  #######################################################
  # Set the weight : later add as option
  lambdas = rep(1, length(dglist))/length(dglist)
  
  #######################################################
  # Compute and Return
  return(median_selection(dglist, lambdas, mymethod, par.iter, par.eps, nCores))
}


# auxiliary functions -----------------------------------------------------
# (0) switching argument
# (1) 2-wasserstein riemannian gradient descent 

# (0) switching argument --------------------------------------------------
#' @keywords internal
#' @noRd
median_selection <- function(dglist, lbd=(rep(1,length(dglist))/length(dglist)), method, 
                                par.iter, par.eps, nCores=0){
  output = switch(method,
                  wass2 = median_wass2(dglist, lbd, par.iter, par.eps, nCores)
                  )
  return(output)
}

# (1) 2-wasserstein riemannian gradient descent ---------------------------
#' @keywords internal
#' @noRd
median_wass2 <- function(dglist, lambdas, par.iter, par.eps, nCores=0){
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
  
  # compute 1 : geometric median for mean vectors
  if (p < 2){
    mout = stats::median(as.vector(mean3))
  } else {
    mean3riem = RiemBase::riemfactory(t(mean3), name="euclidean")
    mout = as.vector(RiemBase::rbase.median(mean3riem)$x)
  }
  
  # compute 2 : covariance
  if (nCores > 1){
    covrun = wass2covs_fmedian_openmp(covs3, lambdas, par.iter, par.eps, nCores)
  } else {
    covrun = wass2covs_fmedian(covs3, lambdas, par.iter, par.eps)
  }
  covout = covrun$x  
  
  # return
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
# # compare with barycenter
# fmean = gauss.barycenter(mylist5d)
# fmed0 = gauss.median(mylist5d)
# fmed3 = gauss.median(mylist5d, nthreads = 3) 
# 
# microbenchmark(
#   fmed0 = gauss.median(mylist5d, nthreads=0),
#   fmed3 = gauss.median(mylist5d, nthreads=3),
#   fmed6 = gauss.median(mylist5d, nthreads=6),
#   fmed9 = gauss.median(mylist5d, nthreads=9)
# )