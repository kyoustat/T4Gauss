#' Pairwise Distance between Gaussian distributions
#' 
#' 
#' @examples 
#' ## test with 1-dimensional Gaussian distributions
#' list1d = list()
#' for (i in 1:20){
#'    val.center  = rnorm(1, sd=0.25)
#'    list1d[[i]] = wrapgauss1d(mean=val.center, sd=1)
#' }
#' for (j in 21:40){
#'    val.center  = rnorm(1, sd=0.25) + 2
#'    list1d[[j]] = wrapgauss1d(mean=val.center, sd=1)
#' }
#' 
#' ## compute pairwise distance 
#' d1wass2 = gauss.pdist(list1d)
#' 
#' ## visualize
#' gcol <- gray((0:64)/64)
#' opar <- par(pty="s")
#' image(d1wass2[,40:1], col=gcol, axes=FALSE, main="wass2 pairwise distance")
#' par(opar)
#' 
#' \dontrun{
#' ## test with 5-dimensional Gaussian distributions
#' list5d = list()
#' for (i in 1:20){
#'    vec.center  = rmvnorm(1, mean=rep(0,5))
#'    list5d[[i]] = wrapgaussNd(mu=vec.center, sigma=diag(5))
#' }
#' mysig = matrix(rnorm(100*5), ncol=5)
#' mysig = t(mysig)%*%mysig
#' for (j in 21:40){
#'    vec.center  = rmvnorm(1, mean=rep(0,5)+1.5)
#'    list5d[[j]] = wrapgaussNd(mu=vec.center, sigma=mysig)
#' }
#' 
#' ## compute pairwise distance
#' d5wass2 = gauss.pdist(list5d)
#' 
#' ## visualize
#' gcol <- gray((0:64)/64)
#' opar <- par(pty="s")
#' image(d5wass2[,40:1], col=gcol, axes=FALSE, main="wass2 pairwise distance")
#' par(opar)
#' }
#' 
#' @export
gauss.pdist <- function(x, y=NULL, method=c("wass2"), as.dist=FALSE){
  #######################################################
  # Preprocessing
  mymethod = match.arg(method)
  mydist   = as.logical(as.dist)
  if (is.list(x)){
    if (!check_list_gauss(x)){
      stop("* gauss.pdist : input 'x' should be a list of 'wrapgauss' objects having same dimension.")
    }
    dglist = x
  } else {
    xcond = base::inherits(x, "wrapgauss")
    ycond = base::inherits(y, "wrapgauss")
    if (!(xcond&&ycond)){
      stop(" gauss.pdist : input 'x' and 'y' should be of class 'wrapgauss'.")
    }
    dglist = list()
    dglist[[1]] = x
    dglist[[2]] = y
    if (!check_list_gauss(dglist)){
      stop("* gauss.pdist : 'x' and 'y' seem to be inconsistent.")
    }
  }
  
  #######################################################
  # Computation with Branching
  output = switch(mymethod,
                  wass2 = gauss.pdist.wass2(dglist)
  )
  
  
  #######################################################
  # Return
  if (mydist){
    return(stats::as.dist(output))
  } else {
    return(output)
  }
}



# auxiliary functions, methods --------------------------------------------
# (1) gauss.pdist.wass2 : 2-wasserstein distance





# (1) gauss.pdist.wass2 ----------------------------------------------------
#' @keywords internal
#' @noRd
gauss.pdist.wass2 <- function(dglist){
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
  
  # compute (don't forget to take a square-root at the end!)
  return(sqrt(wass2_dist(mean3, covs3)))
}