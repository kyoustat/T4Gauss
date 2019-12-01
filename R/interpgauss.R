#' Interpolate Two Gaussian Distributions
#' 
#' method (1-t)N0 + tN1
#' 
#' @examples 
#' ## test with 1-dimensional Gaussian distributions
#' obj1 = wrapgauss1d(mean=1, sd=1)
#' obj5 = wrapgauss1d(mean=5, sd=5)
#' 
#' ## let's interpolate
#' obj.half  = interpgauss(obj1, obj5, t=0.5)
#' 
#' 
#' @export
interpgauss <- function(gauss1, gauss2, t=0.5, method=c("wass2")){
  #######################################################
  # Preprocessing
  if (!inherits(gauss1, "wrapgauss")){
    stop("* interpgauss : an input 'gauss1' should be a S3 object of class 'wrapgauss'.")
  }
  if (!inherits(gauss2, "wrapgauss")){
    stop("* interpgauss : an input 'gauss2' should be a S3 object of class 'wrapgauss'.")
  }
  if ((length(t)>1)||(t>=1)||(t<=0)){
    stop("* interpgauss : 't' should be a real number in (0,1).")
  }
  if (gauss1$dimension!=gauss2$dimension){
    stop("* interpgauss : two 'wrapgauss' objects should be of same dimension.")
  }
  mymethod = match.arg(method)
  
  #######################################################
  # Branching and Execution
  output = switch(mymethod,
                  wass2 = interpgauss_wass2(gauss1, gauss2, t))
  
  #######################################################
  # Return
  return(output)
}


# individual methods ------------------------------------------------------
#' @keywords internal
#' @noRd
interpgauss_wass2 <- function(gauss1, gauss2, t){
  tmpout = wass2_interp(gauss1$mu, gauss2$mu, gauss1$sigma, gauss2$sigma, t)
  mydim  = gauss1$dimension
  if (mydim < 2){
    output = wrapgauss1d(mean=tmpout$mu_t, sd=sqrt(tmpout$sig_t))
  } else {
    output = wrapgaussNd(mu=tmpout$mu_t, sigma=tmpout$sig_t)
  }
  return(output)
}
