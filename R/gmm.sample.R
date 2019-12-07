#' Sample from Gaussian Mixture Model
#' 
#' @export
gmm.sample <- function(n=100, obj.gmm){
  #######################################################
  # Preprocessing
  if (!inherits(obj.gmm, "wrapgmm")){
    stop("* gmm.sample : input 'obj.gmm' should be an object of class 'wrapgmm'.")
  }
  if (obj.gmm$wglist[[1]]$dimension==1){
    mydim = 1
  } else {
    mydim = round(obj.gmm$wglist[[1]]$dimension)
  }
  myweights = obj.gmm$weight
  myclass   = length(obj.gmm$wglist)
  myn       = round(n)
  
  #######################################################
  # sampling index
  sid = base::sample(c(1:myclass, sample(1:myclass, n-myclass, replace=TRUE, prob=myweights)))
  
  #######################################################
  # do the sampling
  #   1. prepare
  if (mydim < 2){
    output = rep(0,myn)
  } else {
    output = array(0,c(myn, mydim))
  }
  #   2. do it
  for (i in 1:myclass){
    sidnow = (sid==i)
    sidn   = sum(sidnow)
    if (mydim < 2){
      output[sidnow]  = gauss.sample(sidn, obj.gmm$wglist[[i]])
    } else {
      output[sidnow,] = gauss.sample(sidn, obj.gmm$wglist[[i]])
    }
  }
  
  #######################################################
  # return
  return(output)
}

# # let's try to make it
# myn   = 50000
# glist = list()
# glist[[1]] = wrapgauss1d(mean=-4, sd=0.5)
# glist[[2]] = wrapgauss1d(mean=0, sd=1)
# glist[[3]] = wrapgauss1d(mean=2, sd=0.5)
# gweight    = rep(1/3, 3)
# gmmobj     = wrapgmm(glist, weight=gweight)
# 
# # 1. sample
# gmmsample  = gmm.sample(n=myn, gmmobj)
# 
# # 2. true density
# draw.x = seq(from=-10,to=10,length.out=1000)
# draw.y = gmm.eval(draw.x, gmmobj)$density
# 
# # visualize
# graphics.off()
# par(mfrow=c(1,2),pty="s")
# plot(draw.x, draw.y, type="l")
# hist(gmmsample, xlim=c(-10,10), freq=TRUE)
