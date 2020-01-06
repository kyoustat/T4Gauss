#' K-Sets Algorithm
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
#' cl2 <- gauss.ksets(mylist, k=2)$cluster
#' cl3 <- gauss.ksets(mylist, k=3)$cluster
#' cl4 <- gauss.ksets(mylist, k=4)$cluster
#' 
#' ## compute 2-dimensional embedding for visualization
#' mds2d <- gauss.mds(mylist, ndim=2)$embed
#' mdsx <- as.vector(mds2d[,1])
#' mdsy <- as.vector(mds2d[,2])
#' 
#' ## visualize
#' opar = par(mfrow=c(1,3), pty="s")
#' plot(mdsx, mdsy, pch=19, col=cl2, main="k=2 means")
#' plot(mdsx, mdsy, pch=19, col=cl3, main="k=3 means")
#' plot(mdsx, mdsy, pch=19, col=cl4, main="k=4 means")
#' par(opar)
#' 
#' \dontrun{
#' ## see the effect of different initialization
#' c3random  = gauss.ksets(mylist, k=3, init.type="random")$cluster
#' c3kpp     = gauss.ksets(mylist, k=3, init.type="kmeans++")$cluster
#' c3medoids = gauss.ksets(mylist, k=3, init.type="kmedoids")$cluster
#' 
#' ## visualize
#' opar = par(mfrow=c(1,3), pty="s")
#' plot(mdsx, mdsy, pch=19, col=c3random,  main="init: random")
#' plot(mdsx, mdsy, pch=19, col=c3kpp,     main="init: k-means++")
#' plot(mdsx, mdsy, pch=19, col=c3medoids, main="init: k-medoids")
#' par(opar)
#' }
#' 
#' @export
gauss.ksets <- function(glist, k=2, type=c("wass2"), maxiter=100, centers=c("mean","median","medoid"),
                        init.type=c("random","kmeans++","kmedoids","specified"), init.label=NULL){
  #######################################################
  # Preprocessing
  if (!check_list_gauss(glist)){
    stop("* gauss.ksets : input 'glist' should be a list of 'wrapgauss' objects having same dimension.")
  }
  mytype = match.arg(type)
  myk    = round(k)
  myn    = length(glist)
  myiter = round(maxiter)
  ctyper = match.arg(centers)

  #######################################################
  # Initialize
  if (myk >= myn){
    stop("* gauss.ksets : 'k' should be a smaller number than 'length(glist)'.")
  }
  init.type = match.arg(init.type)
  if (all(init.type=="random")){
    label.old = base::sample(c(base::sample(1:myk, myn-myk, replace=TRUE), 1:myk))
  } else if (all(init.type=="kmeans++")){
    pdistnow  = gauss.pdist(glist, type=mytype, as.dist=TRUE)
    label.old = DAS::kmeanspp(pdistnow, k=myk)$cluster
  } else if (all(init.type=="kmedoids")){
    label.old = gauss.kmedoids.internal(glist, k=myk)$clustering   
  } else if (all(init.type=="specified")) {
    label.old = as.integer(as.factor(init.label))
    if ((length(label.old)!=myn)||(length(unique(label.old))!=myk)){
      stop("* gauss.ksets : input 'init.label' is not a valid vector of length/number of unique labels.")
    }
  }
  
  #######################################################
  # Plug-in to K-Sets Algorithm
  pdmat = T4Gauss::gauss.pdist(glist, type=mytype)
  outks = DAS::ksets(pdmat, k=myk, init.type="specified", maxiter=myiter, init.label=label.old)
  my.cluster = as.integer(as.factor(outks$cluster))
  ulabel = sort(unique(my.cluster))
  my.centers = list()
  for (i in 1:length(ulabel)){
    idnow = which(my.cluster==ulabel[i])
    if (length(idnow)<2){
      my.centers[[i]] = glist[idnow]
    } else {
      my.centers[[i]] = ksets_selector(glist[idnow], ctyper)
    }
  }
  
  #######################################################
  # Output 
  output = list()
  output$cluster = my.cluster
  output$centers = my.centers
  return(output)
}
  
# selector ----------------------------------------------------------------
#' @keywords internal
#' @noRd
ksets_selector <- function(gdata, ctype){
  if (all(ctype=="mean")){
    return(gauss.barycenter(gdata))
  } else if (all(ctype=="median")){
    return(gauss.median(gdata)) 
  } else if (all(ctype=="medoid")){
    pdmat = gauss.pdist(gdata, type="wass2")
    idmin = which.min(base::rowSums(pdmat))
    return(gdata[idmin])
  }
}