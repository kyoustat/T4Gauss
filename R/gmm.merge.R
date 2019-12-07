#' Large-scale Finite Gaussian Mixture Model Estimation via Split and Merge
#' 
#' 
#' @export
gmm.merge <- function(gmmlist, k=2, ctype=c("kmeans","kmedoids")){
  #######################################################
  # Preprocessing
  if (!check_list_gmm(gmmlist)){
    stop("* gmm.merge : input 'gmmlist' should be a list of 'wrapgauss' objects having same dimension.")
  }
  mynlist = length(gmmlist)
  myk     = round(k)
  myctype = match.arg(ctype)

  #######################################################
  # Rearrange
  #   arr.props (weights)    : mynlist x sum(length(weights))
  #   arr.comps (components) : list of 'wrapgauss' length sum(length(weights))
  for (n in 1:mynlist){
    tgtobj = gmmlist[[n]]
    if (n < 2){
      arr.props = matrix(tgtobj$weight, nrow=1)
      arr.comps = tgtobj$wglist
    } else {
      counter = length(tgtobj$weight)
      arr.comps = c(arr.comps, tgtobj$wglist) # concatenate components
      arr.props = rbind(cbind(arr.props, array(0,c(nrow(arr.props),counter))),c(rep(0,ncol(arr.props)), as.vector(tgtobj$weight)))
    }
  }
  
  #################################################################
  # Step 1. Merge Components via Wasserstein K-Means / K-Medoids
  if (all(myctype=="kmeans")){
    clabel = gauss.kmeans(arr.comps, k=myk)$cluster
  } else {
    clabel = gauss.kmedoids(arr.comps, k=myk)$clustering
  }
  clabel = as.integer(as.factor(clabel))
  clist  = gauss.kmeans.center(arr.comps, clabel, myk, "wass2") # use 'wass2' center
  
  #################################################################
  # Step 2. Merge Component Weights
  ulabel  = unique(clabel)
  cprops  = array(0,c(mynlist, length(ulabel)))
  for (i in 1:length(ulabel)){
    idnow = which(clabel==ulabel[i])
    if (length(idnow)==1){
      cprops[,i] = arr.props[,idnow]
    } else {
      cprops[,i] = base::rowSums(arr.props[,idnow])
    }
  }
  
  #################################################################
  # Step 3. Use RiemSphere to compute mean element of Simplex
  cweight = (as.vector(mle.spnorm(sqrt(cprops), method="Newton")$mu)^2)
  cweight = cweight/sum(cweight)
  
  #################################################################
  # Wrap as GMM object and return
  myobj = wrapgmm(clist, weight=cweight)
  return(myobj)
}

# 
# # test with 'mlbench.smiley'
# library(mlbench)
# library(microbenchmark)
# 
# # generate data
# myk = 5
# myn = 100000
# smiley.large = mlbench.smiley(n=myn)
# smiley.small = list()
# for (i in 1:10){
#   smiley.small[[i]] = mlbench.smiley(n=round(myn/10))
# }
# slx = smiley.large$x
# 
# # run 1. large object
# obj.large = fitgmm(smiley.large$x, k=myk)$gmmobj
# 
# # run 2. small object
# obj.parts = list()
# for (i in 1:10){
#   obj.parts[[i]] = fitgmm(smiley.small[[i]]$x, k=myk)$gmmobj
#   print(paste("iteration ",i,"/10 complete.", sep=""))
# }
# obj.small = gmm.merge(obj.parts, k=myk)
# 
# # let's try visualization
# lab.large = gmm.eval(slx, obj.large)$cluster
# lab.small = gmm.eval(slx, obj.small)$cluster
# 
# par(mfrow=c(1,2))
# plot(slx, pch=19, cex=0.1, col=lab.large, main=paste("fit full for k=",myk,sep=""))
# plot(slx, pch=19, cex=0.1, col=lab.small, main=paste("fit divided for k=",myk,sep=""))



# # 
# x1 = rmvnorm(50, mean=rep(-1,4))
# x2 = rmvnorm(50, mean=rep(+1,4))
# xx = rbind(x1, x2)
# 
# cl1 = fitgmm(xx, k=1)
# cl2 = fitgmm(xx, k=2)
# cl3 = fitgmm(xx, k=3)
# gmmlist = list()
# gmmlist[[1]] = cl1$gmmobj
# gmmlist[[2]] = cl2$gmmobj
# gmmlist[[3]] = cl3$gmmobj
# theobj$wglist = wglist
# theobj$weight = weight
# return(structure(theobj, class="wrapgmm"))