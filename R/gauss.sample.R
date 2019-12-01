#' Draw random samples from Gaussian distribution via 'wrapgauss' object
#' 
#' It's a wrapper for 'rmvnorm' and 'rmvnorm'
#' 
#' @export
gauss.sample <- function(n, wgauss){
  #######################################################
  # Preprocessing
  if (nargs() < 2){
    stop("* gauss.sample : this function requires two arguments, the number of draws and 'wrapgauss' object.")
  }
  if (!inherits(wgauss, "wrapgauss")){
    stop("* gauss.sample : 'wgauss' should be an object of class 'wrapgauss'.")
  }

  #######################################################
  # Case Branching
  mydim = wgauss$dimension
  if (mydim==1){
    return(T4Gauss::rnorm(round(n), mean=mydim$mu, sd=sqrt(mydim$sigma)))
  } else {
    return(T4Gauss::rmvnorm(round(n), mean=mydim$mu, sigma=mydim$sigma))
  }
}