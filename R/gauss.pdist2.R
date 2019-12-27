#' Pairwise Distances between two sets of Gaussian distributions
#' 
#' 
#' @export
gauss.pdist2 <- function(glist1, glist2, type=c("wass2")){
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



# auxiliary functions -----------------------------------------------------
# (0) gauss.pdist2.selector
# (1) gauss.pdist2.wass2 




# (0) gauss.pdist2.selector -----------------------------------------------
#' @keywords internal
#' @noRd
gauss.pdist2.selector <- function(glist1, glist2, mytype){
  output = switch(mytype,
                  wass2 = gauss.pdist2.wass2(glist1, glist2))
  return(output)
}

# (1) gauss.pdist2.wass2 --------------------------------------------------
#' @keywords internal
#' @noRd
gauss.pdist2.wass2 <- function(glist1, glist2){
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
  
  # compute and return (don't forget to sqrt)
  return(sqrt(wass2_dist2(mean1, mean2, covs1, covs2)))
}