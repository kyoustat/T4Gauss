## Auxiliary Function for Computation
#   (1) aux_whichmin     : find the minimal index
#   (2) aux_whichmax     : find the maximal index

# (1) aux_whichmin --------------------------------------------------------
#' @keywords internal
#' @noRd
aux_whichmin <- function(vec){
  mval    = base::min(vec)
  idlarge = which(vec<=mval)
  if (length(idlarge)==1){
    return(idlarge)
  } else {
    return(base::sample(idlarge, 1))
  }
}

# (2) aux_whichmax --------------------------------------------------------
#' @keywords internal
#' @noRd
aux_whichmax <- function(vec){
  mval    = base::max(vec)
  idlarge = which(vec>=mval)
  if (length(idlarge)==1){
    return(idlarge)
  } else {
    return(base::sample(idlarge, 1))
  }
}
