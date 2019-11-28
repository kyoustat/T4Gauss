#' Sampling from Univariate Gaussian Distribution
#' 
#' @export
rnorm <- function(n, mean=0, sd=1){
  #######################################################
  # Preprocessing
  myn    = base::as.integer(n)
  mymean = base::as.double(mean)
  mysd   = base::as.double(sd)
  
  if (!check_number(mymean, pos=FALSE)){
    stop("* rnorm : 'mean' should be a real number.")
  }
  if (!check_number(mysd, pos=TRUE)){
    stop("* rnorm : 'sd' should be a positive real number.")
  }
  
  #######################################################
  # Main Computation
  output = RcppZiggurat::zrnorm(myn)*mysd + mymean
  return(output)
}