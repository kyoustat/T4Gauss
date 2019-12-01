#' Multivariate Gaussian Distribution
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
#' @param x a param for toBar and notToBar
#' @param y a param just for notToBar
#' @return Hard to have one return section for all functions,
#' might want to have a manual list here.
#' @name mvnorm
NULL

#' @rdname mvnorm
#' @export
dmvnorm <- function(x, mean=rep(0,ncol(x)),
                    sigma=diag(length(mean)), log=FALSE){
  #######################################################
  # Preprocessing
  if (is.vector(x)){
    x = base::matrix(x, nrow=1)
  }
  myn = base::as.integer(nrow(x))
  myk = base::ncol(x)
  if (!check_musigma(mean, sigma)){
    stop("* dmvnorm : 'mu' and 'sigma' are not properly valued parameters.")
  }
  if (length(mean)!=myk){
    stop("* dmvnorm : parameters and data have non-matching dimension.")
  }
  mymean  = as.vector(mean)
  mysigma = as.matrix(sigma)
  mylog   = as.logical(log)
  
  #######################################################
  # Main Computation & Return
  if (mylog){
    return(as.vector(cpp_dmvnorm(x, mymean, mysigma)))
  } else {
    return(base::exp(as.vector(cpp_dmvnorm(x, mymean, mysigma))))
  }
}

#' @rdname mvnorm
#' @export
rmvnorm <- function(n, mean=rep(0, nrow(sigma)), 
                    sigma=diag(length(mean)), method=c("eigen","chol")){
  #######################################################
  # Preprocessing
  myn = base::as.integer(n)
  myp = length(mean)
  if (!check_musigma(mean, sigma)){
    stop("* rmvnorm : 'mu' and 'sigma' are not properly valued parameters.")
  }
  mymethod = match.arg(method)
  
  #######################################################
  # Main Computation
  zz = base::matrix(RcppZiggurat::zrnorm(myn*myp), ncol=myp)
  if (all(mymethod=="chol")){
    cholsig = base::chol(sigma)
    scaled  = zz%*%cholsig
  } else {
    eigss  = base::eigen(sigma, symmetric = TRUE)
    lefty  = eigss$vectors %*% diag(as.vector(sqrt(eigss$values)))
    scaled = zz%*%diag(sqrt(eigss$values))%*%t(eigss$vectors)
  }
  output = scaled + base::matrix(rep(mean, myn), ncol=myp, byrow = TRUE)
  
  #######################################################
  # Return
  return(output)
}

# # test : rmvnorm ----------------------------------------------------------
# library(mvtnorm)
# library(microbenchmark)
# 
# myn = 5000
# myp = 10
# mymu  = rnorm(myp)
# A     = matrix(rnorm(10000*myp),ncol=myp)
# mysig = cov(A)
# 
# xe = T4Gauss::rmvnorm(myn, mean=mymu, sigma=mysig, method="eigen")
# xc = T4Gauss::rmvnorm(myn, mean=mymu, sigma=mysig, method="chol")
# xx = mvtnorm::rmvnorm(myn, mean=mymu, sigma=mysig)
# 
# par(mfrow=c(2,2), pty="s")
# image(mysig,   main="original")
# image(cov(xe), main=paste("eigen with d=",round(norm(cov(xe)-mysig),4),sep=""))
# image(cov(xc), main=paste("chol  with d=",round(norm(cov(xc)-mysig),4),sep=""))
# image(cov(xx), main=paste("mvtn  with d=",round(norm(cov(xx)-mysig),4),sep=""))
# 
# stamp <- microbenchmark(
#   xe = T4Gauss::rmvnorm(myn, mean=mymu, sigma=mysig, method="eigen"),
#   xc = T4Gauss::rmvnorm(myn, mean=mymu, sigma=mysig, method="chol"),
#   xx = mvtnorm::rmvnorm(myn, mean=mymu, sigma=mysig)
# )
# graphics.off()
# plot(stamp)
# 
# 
# test : dmvnorm ----------------------------------------------------------
# library(ggplot2)
# library(microbenchmark)
# stamp2 <- microbenchmark(
#   vals1 = T4Gauss::dmvnorm(xe, sigma=mysig, log=FALSE),
#   vals2 = mvtnorm::dmvnorm(xe, sigma=mysig, log=FALSE)
# )
# autoplot(stamp2)



