## Auxiliary Functions
#   (1) check_number     : check univariate number
#   (2) check_musigma    : check mu (mean) and sigma (covariance)
#   (3) check_list_gauss : check whether a list of gaussian distributions
#   (4) check_list_gmm   : check whether a list of GMM models


# (1) check_number --------------------------------------------------------
#' @keywords internal
#' @noRd
check_number <- function(x, pos=TRUE){
  cond1 = ((!is.infinite(x))&&(!is.na(x))&&(length(x)==1))
  if (pos){
    cond2 = base::ifelse(x>0, TRUE, FALSE)
  } else {
    cond2 = TRUE
  }
  if (cond1&&cond2){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# (2) check_musigma -------------------------------------------------------
#' @keywords internal
#' @noRd
check_musigma <- function(x, sigma){
  cond1 = is.vector(x)
  cond2 = (all(!is.infinite(x))&&all(!is.na(x)))
  cond3 = is.matrix(sigma)
  cond4 = (length(x)==nrow(sigma))
  cond5 = base::isSymmetric(sigma, tol=sqrt(.Machine$double.eps))
  if (cond1&&cond2&&cond3&&cond4&&cond5){
    return(TRUE)
  } else {
    if (!cond1){stop("* check_musigma : error 1 : not a vector mean.")}
    if (!cond2){stop("* check_musigma : error 2 : mean contains Inf/NaNs.")}
    if (!cond3){stop("* check_musigma : error 3 : 'sigma' is not a matrix.")}
    if (!cond4){stop("* check_musigma : error 4 : dimensions do not match.")}
    if (!cond5){stop("* check_musigma : error 5 : covariance is not symmetric.")}
    return(FALSE)
  }
}

# (3) check_list_gauss ----------------------------------------------------
#' @keywords internal
#' @noRd
check_list_gauss <- function(wglist){
  extract_dimension <- function(wg){
    return(round(wg$dimension))
  }
  cond1 = is.list(wglist)
  cond2 = (all(unlist(lapply(wglist, inherits, "wrapgauss"))==TRUE))
  cond3 = (length(unique(unlist(lapply(wglist, extract_dimension))))==1)
  if (cond1&&cond2&&cond3){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# (4) check_list_gmm ------------------------------------------------------
#' @keywords internal
#' @noRd
check_list_gmm <- function(gmmlist){
  extract_dimension <- function(myobj){
    return(round(myobj$wglist[[1]]$dimension))
  }
  cond1 = is.list(gmmlist)
  cond2 = (all(unlist(lapply(gmmlist, inherits, "wrapgmm"))==TRUE))
  cond3 = (length(unique(unlist(lapply(gmmlist, extract_dimension))))==1)
  if (cond1&&cond2&&cond3){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

