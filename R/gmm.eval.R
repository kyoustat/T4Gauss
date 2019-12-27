#' Evaluate density and hard clustering of data with GMM model
#' 
#' 
#' @export
gmm.eval <- function(mydata, gmmobj){
  #######################################################
  # Preprocessing
  if (is.vector(mydata)){
    myn   = length(mydata)
    mydim = 1
  } else if (is.matrix(mydata)){
    myn   = nrow(mydata)
    mydim = ncol(mydata)
  }
  if (!inherits(gmmobj, "wrapgmm")){
    stop("* gmm.eval : input 'gmmobj' should be an object of class 'wrapgmm'.")
  }
  if (gmmobj$wglist[[1]]$dimension!=mydim){
    stop("* gmm.eval : dimensions for the data and gmmobj do not match.")
  }
  myweights = gmmobj$weight
  ncomp     = length(myweights)
  
  #######################################################
  # Evaluate
  #   1. density
  mydensity = array(0,c(myn,ncomp))
  for (i in 1:ncomp){
    compi = gmmobj$wglist[[i]]
    if (mydim==1){
      dnow = as.vector(dnorm(mydata, mean=compi$mu, sd=sqrt(compi$sigma), log=FALSE))
    } else {
      dnow = as.vector(dmvnorm(mydata, mean=compi$mu, sigma = compi$sigma, log=FALSE))
    }
    mydensity[,i] = myweights[i]*dnow
  }
  output.density = base::rowSums(mydensity)
  #   2. hard clustering
  output.cluster = apply(mydensity, 1, aux_whichmax)
  
  #######################################################
  # Return
  output = list()
  output$density = output.density
  output$cluster = output.cluster
  return(output)
}