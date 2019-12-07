#' gmmfit
#' 
#' @examples 
#' x1 = rmvnorm(50, mean=rep(-1,4))
#' x2 = rmvnorm(50, mean=rep(+1,4))
#' xx = rbind(x1, x2)
#' 
#' cl1 = fitgmm(xx, k=1)
#' cl2 = fitgmm(xx, k=2)
#' cl3 = fitgmm(xx, k=3)
#' 
#' cls = array(0,c(10,3))
#' for (i in 1:10){
#'   cls[i,] = as.vector(fitgmm(xx, k=i, method="diag")$criteria)
#'   print(paste("iteration ",i," complete..",sep=""))
#' }
#' 
#' yy = c(rnorm(30,mean=-1),rnorm(30,mean=1))
#' cls = array(0,c(10,3))
#' for (i in 1:10){
#'   cls[i,] = as.vector(fitgmm(yy, k=i, method="diag")$criteria)
#'   print(paste("iteration ",i," complete..",sep=""))
#' }
#' 
#' @export
fitgmm <- function(data, k=2, maxiter=10, method=c("full","diag")){ 
  #######################################################
  # Preprocessing
  if (is.vector(data)){
    data = matrix(data, ncol=1)
  }
  mymethod = match.arg(method)
  myk      = round(k)
  myiter   = round(maxiter)
  
  #######################################################
  # Computation
  if (all(mymethod=="full")){
    return(fitgmm.full(data, k=myk, maxiter=myiter))
  } else {
    return(fitgmm.diag(data, k=myk, maxiter=myiter))
  }
}


# case 1. full covariance -------------------------------------------------
#' @keywords internal
#' @noRd
fitgmm.full <- function(data, k=2, maxiter=100){
  myn = nrow(data)
  myp = ncol(data)
  myk = round(k)
  
  tmpout  = arma_gmm_full(t(data), myk, round(maxiter))
  tmpmean = t(tmpout$means)
  tmpcovs = tmpout$covs
  
  # compute and wrap as gmm object
  wglist = list()
  for (k in 1:myk){
    if (myp<2){
      wglist[[k]] = wrapgauss1d(mean=tmpmean[k,], sd=sqrt(tmpcovs[,,k]))
    } else {
      wglist[[k]] = wrapgaussNd(mu=tmpmean[k,], sigma=tmpcovs[,,k])
    }
  }
  myweight = as.vector(tmpout$weight)
  myweight = myweight/sum(myweight)
  gmmobj   = wrapgmm(wglist, weight=myweight)

  # score part
  loglkd = tmpout$loglkd
  par.k  = (myp*myk) + ((myp*(myp+1)/2)*myk) + (myk-1) # mean + covs + proportion  
  
  gmmscore = rep(0,3)
  gmmscore[1] = -2*loglkd + (2*par.k)             # AIC  : minimum
  gmmscore[2] = -2*loglkd + (par.k*log(myn))      # BIC  : minimum
  gmmscore[3] = -2*loglkd + 2*par.k*log(log(myn)) # HQIC : minimum
  gmmscore    = matrix(gmmscore, nrow=1)
  colnames(gmmscore) = c("AIC","BIC","HQIC")
  
  # return
  output = list()
  output$gmmobj = gmmobj
  output$criteria = gmmscore
  return(output)
}


# case 2. diag covariance -------------------------------------------------
#' @keywords internal
#' @noRd
fitgmm.diag <- function(data, k=2, maxiter=100){
  myn = nrow(data)
  myp = ncol(data)
  myk = round(k)
  
  tmpout  = arma_gmm_diag(t(data), myk, round(maxiter))
  tmpmean = t(tmpout$means) # rows are means
  tmpcovs = t(tmpout$covs)  # rows are diagonals
  
  # compute and wrap as gmm object
  wglist = list()
  for (k in 1:myk){
    if (myp<2){
      wglist[[k]] = wrapgauss1d(mean=tmpmean[k,], sd=sqrt(tmpcovs[k,]))
    } else {
      wglist[[k]] = wrapgaussNd(mu=tmpmean[k,], sigma=diag(tmpcovs[k,]))
    }
  }
  myweight = as.vector(tmpout$weight)
  myweight = myweight/sum(myweight)
  gmmobj   = wrapgmm(wglist, weight=myweight)
  
  # score part
  loglkd = tmpout$loglkd
  par.k  = (myp*myk) + (myp*myk) + (myk-1) # mean + covs + proportion  
  
  gmmscore = rep(0,3)
  gmmscore[1] = -2*loglkd + (2*par.k)             # AIC  : minimum
  gmmscore[2] = -2*loglkd + (par.k*log(myn))      # BIC  : minimum
  gmmscore[3] = -2*loglkd + 2*par.k*log(log(myn)) # HQIC : minimum
  gmmscore    = matrix(gmmscore, nrow=1)
  colnames(gmmscore) = c("AIC","BIC","HQIC")
  
  # return
  output = list()
  output$gmmobj = gmmobj
  output$criteria = gmmscore
  return(output)
}

