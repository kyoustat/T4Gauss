#' MDS for Gaussian
#' 
#' 
#' @export
gauss.mds <- function(glist, ndim=2, type=c("wass2")){
  #######################################################
  # Preprocessing
  if (!check_list_gauss(glist)){
    stop("* gauss.mds : input 'glist' should be a list of 'wrapgauss' objects having same dimension.")
  }
  mydistance = match.arg(type)
  myndim     = round(ndim)
  
  #######################################################
  # Compute Pairwise Distance
  pdmat = stats::as.dist(gauss.pdist.selector(glist, mytype=mydistance))
  output = DAS::cmds(pdmat, ndim=myndim)
  
  #######################################################
  # Return
  return(output)
}