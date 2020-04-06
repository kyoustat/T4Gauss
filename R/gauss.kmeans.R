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
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(3,1))
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
                         init.type=c("random","kmeans++","kmedoids","specified"), init.label=NULL,
                         permute.order=FALSE){
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
  if (missing(init.type)){
    init.type = "random"
  } else {
    init.type = match.arg(init.type)  
  }
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
  # center.old = gauss.kmeans.center(glist, label.old, myk, mytype, nCores)
  
  ############################################################
  # Re-coding MacQueen-type Algorithm
  # 1. renaming some parameters
  K = myk
  N = myn
  maxiter = round(maxiter)
  
  # 2. label information
  label = label.old
  index = list()
  for (k in 1:K){
    index[[k]] = which(label==k)
  }
  
  # 3. Initialization
  # 3-1. centroids
  list.centroids = list()
  for (k in 1:K){
    list.centroids[[k]] = gkmeans.center(glist, index[[k]], mytype, nCores)
  }
  # 3-2. SSE 
  SSEold = gkmeans.SSE(glist, list.centroids, index, mytype)
  print(paste0("iteration 0"," & SSE=",SSEold))
  
  # 4. Main Iteration
  for (it in 1:maxiter){
    # 4-1. update for each element 
    if (permute.order){
      vecn = base::sample(N)
    } else {
      vecn = 1:N
    }
    for (n in vecn){ # 1:N or random update
      if (length(index[[label[n]]]) > 1){ # Idea : If a Singleton Set, don't need to update
        old.class = label[n]
        now.xlist = glist[n]
        now.d2c   = gkmeans.D2Centroids(now.xlist, list.centroids, mytype)
        now.class = which.min(now.d2c)
        
        if (now.class!=old.class){ 
          label[n] = now.class
          index[[old.class]] = setdiff(index[[old.class]], n)
          index[[now.class]] = c(index[[now.class]], n)
          
          list.centroids[[old.class]] = gkmeans.center(glist, index[[old.class]], mytype, nCores)
          list.centroids[[now.class]] = gkmeans.center(glist, index[[now.class]], mytype, nCores) 
        } 
      }
    }
    # 4-2. compute SSE
    SSEnew = gkmeans.SSE(glist, list.centroids, index, mytype)
    print(paste0("iteration ",it," & SSE=",SSEnew))

    # 4-3. update
    SSEinc = abs(SSEold-SSEnew)/SSEold
    SSEold = SSEnew
    if ((SSEinc < 1e-5)&&(it > 2)){
      break
    }
  }

  ############################################################
  # Return
  outlab = rep(0,N)
  for (k in 1:K){
    outlab[index[[k]]] = k
  }

  output = list()
  output$cluster = outlab
  output$means   = list.centroids
  return(output) 
}


# new set of auxiliary functions ------------------------------------------
#' @export
gkmeans.center <- function(glist, indexvec, type, nthreads){
  if (length(indexvec)==1){
    output = glist[[indexvec]]
  } else {
    gparts = glist[indexvec]
    if (nthreads > 1){
      output = barygauss_selection(gparts, method="wass2rgd", 
                                   par.iter=100, par.eps=1e-6, nCores=nthreads)  
    } else {
      output = barygauss_selection(gparts, method="wass2fpt", 
                                   par.iter=100, par.eps=1e-6) 
    }
  }
  return(output)
}
#' @export
gkmeans.SSE <- function(glist, list.centroids, index, type){
  K   = length(index)
  SSE = 0
  for (k in 1:K){
    dvec = as.vector(gauss.pdist2.selector(list.centroids[k], glist[index[[k]]], type))
    SSE  = SSE + sum(dvec^2)
  }
  return(SSE)
}
#' @export
gkmeans.D2Centroids <- function(list.xnow, list.centroids, type){
  return(as.vector(gauss.pdist2.selector(list.xnow, list.centroids, type)))
}