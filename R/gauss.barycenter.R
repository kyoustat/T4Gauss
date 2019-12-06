#' Barycenter of Gaussian Distributions
#' 
#' @export
gauss.barycenter <- function(x, y=NULL, type=c("wass2")){
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
  
  
  #######################################################
  # Set the weight
  lambdas = rep(1, length(dglist))/length(dglist)
  
  #######################################################
  # Compute and Return
  return(barygauss_selection(dglist, lambdas, mymethod))
}



# auxiliary functions -----------------------------------------------------
# (0) switching argument
# (1) type : wass2

# (0) switching argument --------------------------------------------------
#' @keywords internal
#' @noRd
barygauss_selection <- function(dglist, lbd=(rep(1,length(dglist))/length(dglist)), method){
  output = switch(method,
                  wass2 = barygauss_wass2(dglist, lbd))
  return(output)
}

# (1) wass2 ---------------------------------------------------------------
#' @keywords internal
#' @noRd
barygauss_wass2 <- function(dglist, lambdas){
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
  tmpout = wass2_barycenter(mean3, covs3, lambdas)
  if (p < 2){
    return(wrapgauss1d(mean=tmpout$mu_bary, sd=sqrt(tmpout$sig_bary)))
  } else {
    return(wrapgaussNd(mu=tmpout$mu_bary, sigma=tmpout$sig_bary))
  }
}