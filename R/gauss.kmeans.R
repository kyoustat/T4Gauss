#' K-Means clustering for Gaussian Distributions
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
#' cl2 <- gauss.kmeans(mylist, k=2)$cluster
#' cl3 <- gauss.kmeans(mylist, k=3)$cluster
#' cl4 <- gauss.kmeans(mylist, k=4)$cluster
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
#' c3random  = gauss.kmeans(mylist, k=3, init.type="random")$cluster
#' c3kpp     = gauss.kmeans(mylist, k=3, init.type="kmeans++")$cluster
#' c3medoids = gauss.kmeans(mylist, k=3, init.type="kmedoids")$cluster
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
gauss.kmeans <- function(glist, k=2, type=c("wass2"), maxiter=100, nthreads=1,
                         init.type=c("random","kmeans++","kmedoids","specified"), init.label=NULL){
  #######################################################
  # Preprocessing
  if (!check_list_gauss(glist)){
    stop("* gauss.kmeans : input 'glist' should be a list of 'wrapgauss' objects having same dimension.")
  }
  mytype = match.arg(type)
  myk    = round(k)
  myn    = length(glist)
  nCores = round(nthreads)
  
  #######################################################
  # Initialize
  if (myk >= myn){
    stop("* gauss.kmeans : 'k' should be a smaller number than 'length(glist)'.")
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
      stop("* gauss.kmeans : input 'init.label' is not a valid vector of length/number of unique labels.")
    }
  }
  center.old = gauss.kmeans.center(glist, label.old, myk, mytype, nCores)
  
  #######################################################
  # Naive Algorithm
  for (it in 1:maxiter){
    # Assignment Step
    # A-1. compute pairwise distance (N x K)
    pdmat = gauss.pdist2(glist, center.old, type=mytype)
    # A-2. class assignment
    label.new = as.integer(as.factor(base::apply(pdmat, 1, aux_whichmin)))
    label.new = gauss.kmeans.label.adjust(glist, label.new, myk)
  
    # Update Step
    center.new = gauss.kmeans.center(glist, label.new, myk, mytype, nCores)
  
    # Iteration Control
    labeldel   = as.double(mclustcomp::mclustcomp(label.new, label.old,types="nmi1")[2])
    label.old  = label.new
    center.old = center.new
    if ((labeldel>=0.99)&&(it>=5)){
      break
    }
  }
  
  ############################################################
  # Return
  output = list()
  output$cluster = label.old
  output$means   = center.old
  return(output)
}

  

# auxiliary functions -----------------------------------------------------
#' @keywords internal
#' @noRd
gauss.kmeans.center <- function(glist, label, k, type, nthreads){
  label  = round(label)
  n      = length(glist)
  k      = round(k)

  centers = list()
  for (i in 1:k){
    idnow = which(label==i)
    if (length(idnow)==1){
      centers[[i]] = glist[[idnow]]
    } else {
      gparts = glist[idnow]
      if (all(type=="wass2")){
        if (nthreads > 1){
          centers[[i]] = barygauss_selection(gparts, method="wass2rgd", 
                                             par.iter=100, par.eps=1e-6, nCores=nthreads)  
        } else {
          centers[[i]] = barygauss_selection(gparts, method="wass2fpt", 
                                             par.iter=100, par.eps=1e-6) 
        }
      }
    }
  }
  return(centers)
}

#' @keywords internal
#' @noRd
gauss.kmeans.label.adjust <- function(glist, label, k){ # if
  rk = round(k)
  lul = length(unique(label))
  if (lul==rk){
    return(label)
  } else {
    newlab = label
    ids    = which(label == as.integer(names(which.max(table(label)))))
    gpart  = glist[ids]
    newlab[ids] = gauss.kmedoids.internal(gpart, k=(rk-lul+1))$clustering + rk
    return(as.integer(as.factor(newlab)))
  }
}