#' Draw random samples from Gaussian distribution via 'wrapgauss' object
#' 
#' It's a wrapper for 'rmvnorm' and 'rmvnorm'
#' 
#' @export
gauss.sample <- function(n=100, obj.gauss){
  #######################################################
  # Preprocessing
  if (nargs() < 2){
    stop("* gauss.sample : this function requires two arguments, the number of draws and 'wrapgauss' object.")
  }
  if (!inherits(obj.gauss, "wrapgauss")){
    stop("* gauss.sample : 'obj.gauss' should be an object of class 'wrapgauss'.")
  }

  #######################################################
  # Case Branching
  mydim = obj.gauss$dimension
  if (mydim==1){
    return(T4Gauss::rnorm(round(n), mean=obj.gauss$mu, sd=sqrt(obj.gauss$sigma)))
  } else {
    return(T4Gauss::rmvnorm(round(n), mean=obj.gauss$mu, sigma=obj.gauss$sigma))
  }
}