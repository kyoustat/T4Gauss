#' kmedoids
#' 
#' k-Medoids is a generally applicable clustering algorithm 
#' as long as we have concept of dissimilarity. We adopt \code{pam} algorithm 
#' by \pkg{cluster} package. See \code{\link[cluster]{pam}} for more details.
#' 
#' @param glist list of objects, a S3 object of \code{riemdata} class. See \code{\link{riemfactory}} for more details.
#' @param k the number of clusters.
#' @param type type of distance metric to be used.
#' 
#' 
#' @examples 
#' ## generate three-cluster data with univariate Gaussians
#' mylist = list()
#' for (i in 1:10){
#'    mylist[[i]] = wrapgauss1d(mean=-2-runif(1), sd=runif(1))
#' }
#' for (i in 11:20){
#'    mylist[[i]] = wrapgauss1d(mean=0, sd=runif(1))
#' }
#' for (i in 21:30){
#'    mylist[[i]] = wrapgauss1d(mean=2+runif(1), sd=runif(1))
#' }
#' 
#' ## apply clustering with different k values
#' cl2 <- gauss.kmedoids(mylist, k=2)$cluster
#' cl3 <- gauss.kmedoids(mylist, k=3)$cluster
#' cl4 <- gauss.kmedoids(mylist, k=4)$cluster
#' 
#' ## compute 2-dimensional embedding for visualization
#' mds2d <- gauss.mds(mylist, ndim=2)$embed
#' mdsx <- as.vector(mds2d[,1])
#' mdsy <- as.vector(mds2d[,2])
#' 
#' ## visualize
#' opar = par(mfrow=c(1,3), pty="s")
#' plot(mdsx, mdsy, pch=19, col=cl2, main="k=2 medoids")
#' plot(mdsx, mdsy, pch=19, col=cl3, main="k=3 medoids")
#' plot(mdsx, mdsy, pch=19, col=cl4, main="k=4 medoids")
#' par(opar)
#' 
#' @export
gauss.kmedoids <- function(glist, k=2, type=c("wass2")){
  #######################################################
  # Preprocessing
  if (!check_list_gauss(glist)){
    stop("* gauss.kmedoids : input 'glist' should be a list of 'wrapgauss' objects having same dimension.")
  }
  mytype = match.arg(type)
  myk    = round(k)
  
  #######################################################
  # Compute Pairwise Distance and run PAM
  tmpdist = gauss.pdist(glist, type=mytype, as.dist=TRUE)
  pamout  = cluster::pam(tmpdist, k=myk)
  
  #######################################################
  # Return
  output = list()
  output$cluster = pamout$clustering
  output$medoids = glist[pamout$id.med]
  return(output)
}



# internal function : default as 2-wasserstein distance -------------------
#' @keywords internal
#' @noRd
gauss.kmedoids.internal <- function(glist, k=2){
  #######################################################
  # Preprocessing
  myk    = round(k)
  
  #######################################################
  # Compute Pairwise Distance and run PAM
  tmpdist = gauss.pdist(glist, type="wass2", as.dist=TRUE)
  pamout  = cluster::pam(tmpdist, k=myk)
  
  #######################################################
  # Return
  return(pamout)
}

