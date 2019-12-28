#' Pairwise Distance/Dissimilarities between Gaussian distributions
#' 
#' 
#' @examples 
#' ## test with 1-dimensional Gaussian distributions
#' list1d = list()
#' for (i in 1:15){
#'    val.center  = rnorm(1, sd=0.25)
#'    val.sd      = runif(1, min=0.75,max=1)
#'    list1d[[i]] = wrapgauss1d(mean=val.center, sd=val.sd)
#' }
#' for (j in 16:40){
#'    val.center  = rnorm(1, sd=0.25) + 10
#'    val.sd      = runif(1, min=4.75,max=5)
#'    list1d[[j]] = wrapgauss1d(mean=val.center, sd=val.sd)
#' }
#' 
#' ## compute pairwise distance 
#' d1wass2 = gauss.pdist(list1d, type="wass2")
#' d1kl    = gauss.pdist(list1d, type="kl")
#' d1skl   = gauss.pdist(list1d, type="skl")
#' d1cs    = gauss.pdist(list1d, type="cs")
#' d1bh    = gauss.pdist(list1d, type="bh")
#' 
#' ## visualize
#' opar <- par(mfrow=c(2,3), pty="s")
#' image(d1wass2[,40:1], axes=FALSE, main="pdist: 2-Wasserstein")
#' image(d1kl[,40:1],    axes=FALSE, main="pdist: KL")
#' image(d1skl[,40:1],   axes=FALSE, main="pdist: Symmetric KL")
#' image(d1cs[,40:1],    axes=FALSE, main="pdist: Cauchy-Schwarz")
#' image(d1bh[,40:1],    axes=FALSE, main="pdist: Bhattacharyya")
#' par(opar)
#' 
#' \dontrun{
#' ## test with 5-dimensional Gaussian distributions
#' list5d = list()
#' for (i in 1:20){
#'    vec.center  = rmvnorm(1, mean=rep(0,5))
#'    list5d[[i]] = wrapgaussNd(mu=vec.center, sigma=diag(5))
#' }
#' for (j in 21:40){
#'    vec.center  = rmvnorm(1, mean=rep(0,5)+1.5)
#'    mysig = matrix(rnorm(100*5), ncol=5)
#'    mysig = t(mysig)%*%mysig
#'    list5d[[j]] = wrapgaussNd(mu=vec.center, sigma=mysig)
#' }
#' 
#' ## compute pairwise distance
#' d5wass2 = gauss.pdist(list5d, type="wass2")
#' d5kl    = gauss.pdist(list5d, type="kl")
#' d5skl   = gauss.pdist(list5d, type="skl")
#' d5cs    = gauss.pdist(list5d, type="cs")
#' d5bh    = gauss.pdist(list5d, type="bh")
#' 
#' ## visualize
#' opar <- par(mfrow=c(2,3), pty="s")
#' image(d5wass2[,40:1], axes=FALSE, main="pdist: 2-Wasserstein")
#' image(d5kl[,40:1],    axes=FALSE, main="pdist: KL")
#' image(d5skl[,40:1],   axes=FALSE, main="pdist: Symmetric KL")
#' image(d5cs[,40:1],    axes=FALSE, main="pdist: Cauchy-Schwarz")
#' image(d5bh[,40:1],    axes=FALSE, main="pdist: Bhattacharyya")
#' par(opar)
#' }
#' 
#' @export
gauss.pdist <- function(x, y=NULL, type=c("bh","cs","kl","skl","wass2"), as.dist=FALSE){
  #######################################################
  # Preprocessing
  mymethod = match.arg(type)
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
  output = gauss.pdist.selector(dglist, mymethod)
  
  #######################################################
  # Return
  if (mydist){
    return(stats::as.dist(output))
  } else {
    return(output)
  }
}





# Auxiliary : selector ----------------------------------------------------
#' @keywords internal
#' @noRd
gauss.pdist.selector <- function(dglist, mytype){
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
  
  # switching
  if (all(mytype=="skl")){
    tmpout = kl_dist(mean3, covs3)
    output = tmpout + t(tmpout)
  } else {
    output = switch(mytype,
                    wass2 = sqrt(wass2_dist(mean3, covs3)),
                    kl    = kl_dist(mean3, covs3),
                    cs    = cs_dist(mean3, covs3),
                    bh    = bh_dist(mean3, covs3))
  }

  return(output)
}

# Personal Notes ----------------------------------------------------------
# wass2    : from OMT Theory of GMM
# kl & skl : Spurek.2016 : Kullbeck-Leibler Divergence (skl = kl(p,q) + kl(q,p))
# cs       : Spurek.2016 : Cauchy-Schwarz Divergence
# bh       : Spurek.2016 : Bhattacharyya Distance
