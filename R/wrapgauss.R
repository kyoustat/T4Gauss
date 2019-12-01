#' Wrap Gaussian distribution into \code{S3} object of class 'wrapgauss'
#' 
#' Group of functions Description section
#' 
#' Group of functions Details paragraph.
#'
#' @section After Arguments and Value sections:
#' Despite its location, this actually comes after the Arguments and Value sections.
#' Also, don't need to use null, could annotate first function, and then
#' using function name as the groupBy name is more intuitive.
#' 
#' @param mean 1d
#' @param sd 1d
#' @param mu nd
#' @param sigma nd
#' 
#' 
#' 
#' @name wrapgauss
NULL

#' @rdname wrapgauss
#' @export
wrapgauss1d <- function(mean, sd){
  #######################################################
  # Preprocessing
  if (nargs() < 2){
    stop("* wrapgauss1d : this function needs exactly two variables 'mean' and 'sd'.")
  }
  mymean = base::as.double(mean)
  mysd   = base::as.double(sd)
  
  if (!check_number(mymean, pos=FALSE)){
    stop("* wrapgauss1d : 'mean' should be a real number.")
  }
  if (!check_number(mysd, pos=TRUE)){
    stop("* wrapgauss1d : 'sd' should be a positive real number.")
  }
  
  heymean  = as.vector(mymean)
  heysigma = base::matrix(mysd^2, ncol=1)
  
  #######################################################
  # Wrap & Return
  output = list()
  output$mu    = heymean
  output$sigma = heysigma
  output$dimension = 1
  
  return(structure(output, class="wrapgauss"))
}

#' @rdname wrapgauss
#' @export
wrapgaussNd <- function(mu, sigma){
  #######################################################
  # Preprocessing
  if (nargs() < 2){
    stop("* wrapgaussNd : this function needs exactly two variables 'mu' and 'sigma'.")
  }
  mymean  = as.vector(mu)
  mysigma = as.matrix(sigma)
  if (!check_musigma(mymean, mysigma)){
    stop("* wrapgaussNd : 'mu' and 'sigma' are not properly valued parameters.")
  }
  if (min(base::eigen(mysigma, symmetric = TRUE, only.values = TRUE)$values) <= 0){
    stop("* wrapgaussNd : 'sigma' should be a positive definite matrix.")
  }
  
  #######################################################
  # Wrap & Return
  output = list()
  output$mu    = mymean
  output$sigma = mysigma
  output$dimension = nrow(mysigma)
  
  return(structure(output, class="wrapgauss"))
}









