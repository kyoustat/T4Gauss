#' Interpolate Two Gaussian Distributions
#' 
#' method (1-t)N0 + tN1
#' 
#' @examples 
#' ## test with 1-dimensional Gaussian distributions
#' objA = wrapgauss1d(mean=1, sd=0.5)
#' objZ = wrapgauss1d(mean=5, sd=2)
#' 
#' ## let's interpolate for t=0.1 to t=0.9
#' vec.t   = seq(from=0.1, to=0.9, by=0.1)
#' vec.obj = list()
#' for (i in 1:9){
#'    vec.obj[[i]] = gauss.interp(objA, objZ, t=vec.t[i])
#' }
#' 
#' ## let's visualize
#' draw.x = seq(from=-3, to=11, length.out=1234)
#' draw.A = gauss.eval(draw.x, objA)
#' draw.Z = gauss.eval(draw.x, objZ)
#' 
#' plot(draw.x, draw.A, type="l", xlim=c(-1,11), ylim=c(0,1), lwd=2,
#'      xlab="x", ylab="y", main="Interpolated Gaussians between Two Black")
#' lines(draw.x, draw.Z, lwd=2)
#' for (i in 1:9){
#'   draw.t = gauss.eval(draw.x, vec.obj[[i]])
#'   lines(draw.x, draw.t, col=(i+1), lty=2, lwd=0.5)
#' }
#' 
#' 
#' @export
gauss.interp <- function(gauss1, gauss2, t=0.5, method=c("wass2")){
  #######################################################
  # Preprocessing
  if (!inherits(gauss1, "wrapgauss")){
    stop("* gauss.interp : an input 'gauss1' should be a S3 object of class 'wrapgauss'.")
  }
  if (!inherits(gauss2, "wrapgauss")){
    stop("* gauss.interp : an input 'gauss2' should be a S3 object of class 'wrapgauss'.")
  }
  if ((length(t)>1)||(t>=1)||(t<=0)){
    stop("* gauss.interp : 't' should be a real number in (0,1).")
  }
  if (gauss1$dimension!=gauss2$dimension){
    stop("* gauss.interp : two 'wrapgauss' objects should be of same dimension.")
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
