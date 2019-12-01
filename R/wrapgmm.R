#' Wrap Mixture of Gaussians into \code{S3} object of class 'wrapgmm'
#' 
#' @export
wrapgmm <- function(wglist, weight=rep(1/length(wglist), length(wglist))){
  #######################################################
  # Preprocessing
  #   1. should be a list
  if (!is.list(wglist)){
    stop("* wrapgmm : 'wglist' should be a list.")
  }
  #   2. all from 'wrapgauss' object
  if (!(all(unlist(lapply(wglist, inherits, "wrapgauss"))==TRUE))){
    stop("* wrapgmm : 'wglist' should be a list of 'wrapgauss' objects.")
  }
  #   3. all have same dimension
  extract_dimension <- function(wg){
    return(round(wg$dimension))
  }
  if (length(unique(unlist(lapply(wglist, extract_dimension))))!=1){
    stop("* wrapgmm : 'wglist' should be a list of 'wrapgauss' objects having same dimension.")
  }
  #   4. weight
  if ((!is.vector(weight))||(sum(weight)!=1)||(length(weight)!=length(wglist))){
    stop("* wrapgmm : 'weight' should be a vector of same length as 'wglist' which sums to 1.")
  }
  
  #######################################################
  # Wrap it
  theobj = list()
  theobj$wglist = wglist
  theobj$weight = weight
  return(structure(theobj, class="wrapgmm"))
}