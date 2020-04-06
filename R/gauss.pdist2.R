#' Pairwise Distances between two sets of Gaussian distributions
#' 
#' 
#' @examples 
#' ## test with 1-dimensional Gaussian distributions
#' list1dA = list()
#' list1dB = list()
#' for (i in 1:5){
#'    val.center  = rnorm(1, sd=0.25)
#'    val.sd      = runif(1, min=0.75,max=1)
#'    list1dA[[i]] = wrapgauss1d(mean=val.center, sd=val.sd)
#' }
#' for (i in 6:10){
#'    val.center  = rnorm(1, sd=0.25) + 10
#'    val.sd      = runif(1, min=4.75,max=5)
#'    list1dA[[i]] = wrapgauss1d(mean=val.center, sd=val.sd)
#' }
#' for (j in 1:11){
#'    val.center   = rnorm(1, sd=0.25)
#'    val.sd       = runif(1, min=0.75,max=1)
#'    list1dB[[j]] = wrapgauss1d(mean=val.center, sd=val.sd)
#' }
#' for (j in 12:25){
#'    val.center  = rnorm(1, sd=0.25) + 10
#'    val.sd      = runif(1, min=4.75,max=5)
#'    list1dB[[j]] = wrapgauss1d(mean=val.center, sd=val.sd)
#' }
#' 
#' ## compute pairwise distance 
#' d1wass2 = gauss.pdist2(list1dA, list1dB, type="wass2")
#' d1kl    = gauss.pdist2(list1dA, list1dB, type="kl")
#' d1skl   = gauss.pdist2(list1dA, list1dB, type="skl")
#' d1cs    = gauss.pdist2(list1dA, list1dB, type="cs")
#' d1bh    = gauss.pdist2(list1dA, list1dB, type="bh")
#' 
#' ## visualize
#' opar <- par(mfrow=c(2,3))
#' image(d1wass2, axes=FALSE, main="pdist2: 2-Wasserstein")
#' image(d1kl,    axes=FALSE, main="pdist2: KL")
#' image(d1skl,   axes=FALSE, main="pdist2: Symmetric KL")
#' image(d1cs,    axes=FALSE, main="pdist2: Cauchy-Schwarz")
#' image(d1bh,    axes=FALSE, main="pdist2: Bhattacharyya")
#' par(opar)
#' 
#' @export
gauss.pdist2 <- function(glist1, glist2, type=c("bh","cs","kl","skl","wass2")){
  #######################################################
  # Preprocessing
  mymethod = match.arg(type)
  xcond = all(unlist(lapply(glist1, inherits, "wrapgauss"))==TRUE)
  ycond = all(unlist(lapply(glist2, inherits, "wrapgauss"))==TRUE)
  if (!(xcond&&ycond)){
    stop(" gauss.pdist2 : input 'glist1' and 'glist2' should be of class 'wrapgauss'.")
  }
  if (!check_list_gauss(glist1)){
    stop("* gauss.pdist2 : 'glist1' is not a valid list of Gaussian distributions of same dimension.")
  }
  if (!check_list_gauss(glist2)){
    stop("* gauss.pdist2 : 'glist2' is not a valid list of Gaussian distributions of same dimension.")
  }
  if (glist1[[1]]$dimension!=glist2[[1]]$dimension){
    stop("* gauss.pdist2 : 'glist1' and 'glist2' are not of same dimension.")
  }
  
  #######################################################
  # Compute and Return
  return(gauss.pdist2.selector(glist1, glist2, mymethod))
}




# Auxiliary : Selector ----------------------------------------------------
#  List of Gaussian Objects
#' @keywords internal
#' @noRd
gauss.pdist2.selector <- function(glist1, glist2, mytype){
  # parameters
  M = length(glist1)
  N = length(glist2)
  p = glist1[[1]]$dimension
  
  # stack for glist1
  mean1 = array(0,c(M,p))
  covs1 = array(0,c(p,p,M))
  for (m in 1:M){
    mean1[m,]  = as.vector(glist1[[m]]$mu)
    covs1[,,m] = as.matrix(glist1[[m]]$sigma)
  }
  
  # stack for glist2
  mean2 = array(0,c(N,p))
  covs2 = array(0,c(p,p,N))
  for (n in 1:N){
    mean2[n,]  = as.vector(glist2[[n]]$mu)
    covs2[,,n] = as.matrix(glist2[[n]]$sigma)
  }
  
  # compute and return
  output = switch(mytype,
                  wass2 = sqrt(wass2_dist2(mean1, mean2, covs1, covs2)),
                  kl    = kl_dist2(mean1, covs1, mean2, covs2),
                  skl   = skl_dist2(mean1, covs1, mean2, covs2),
                  cs    = cs_dist2(mean1, covs1, mean2, covs2),
                  bh    = bh_dist2(mean1, covs1, mean2, covs2))
  return(output)
}
