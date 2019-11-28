## Auxiliary Functions
#   (1) check_number  : check univariate number
#   (2) check_musigma : check mu (mean) and sigma (covariance)



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
  cond5 = base::isSymmetric(sigma)
  if (cond1&&cond2&&cond3&&cond4&&cond5){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
