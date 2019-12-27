#' Evaluate Gaussian
#' 
#' @export
gauss.eval <- function(mydata, obj.gauss, log=FALSE){
  #######################################################
  # Preprocessing
  name.obj  = base::deparse(base::substitute(obj.gauss))
  name.data = base::deparse(base::substitute(mydata))
  
  if (is.vector(mydata)){
    myn   = length(mydata)
    mydim = 1
  } else if (is.matrix(mydata)){
    myn   = nrow(mydata)
    mydim = ncol(mydata)
  }
  if (!inherits(obj.gauss, "wrapgauss")){
    stop(paste("* gauss.eval : input '",name.obj,"' should be an object of class 'wrapgauss'."))
  }
  if (obj.gauss$dimension!=mydim){
    stop(paste("* gauss.eval :'",name.data,"' and '",name.obj,"' are of different dimension."))
  }
  mylog = as.logical(log)
  
  #######################################################
  # do it
  if (mydim < 2){
    dnow = as.vector(dnorm(mydata, mean=obj.gauss$mu, sd=sqrt(obj.gauss$sigma), log=mylog))
  } else {
    dnow = as.vector(dmvnorm(mydata, mean=obj.gauss$mu, sigma = obj.gauss$sigma, log=mylog))
  }
  
  #######################################################
  # return
  return(dnow)
}
