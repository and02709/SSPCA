#' This function decomposes a matrix using singular value decomposition.  
#' The resulting eigenvector corresponding to the loadings in then shrunk
#' using an L1 penalty.
#' @param x this is the matrix which is decomposed and the decomposition vectors
#'   corresponding to the loadings are shrunk
#' @param sumabsu this parameter determines by how much the decomposition vectors
#'   corresponding to the observations
#' @param sumabsv this parameter determines by how much the decomposition vectors
#'   corresponding to the loadings
#' @param niter number of iterations to allow for convergence
#' @param trace option to deisplay potential warnings
#' @param v allows a particular loadings vector to be used
#' @param upos allows the option to force the vector corresponding to 
#'   observations to correspond to the parallel maxima using pmax
#' @param uneg allows the option to force the vector corresponding to 
#'   observations to correspond to the parallel minima using pmin
#' @param vpos allows the option to force the vector corresponding to 
#'   loadings to correspond to the parallel maxima using pmax
#' @param vneg allows the option to force the vector corresponding to 
#'   loadings to correspond to the parallel minima using pmin
#' @keywords singular matrix
#' @export
#' @examples SMD(x, sumabsu, sumabsv, niter=20,trace=TRUE, v, upos, uneg, vpos, vneg)



# Function accepts the matrix x to decompose and employs sparseness using
#   L1 penalty on the vectors corresponding the loadings
SMD <- function(x, sumabsu, sumabsv, niter=20,trace=TRUE, v, upos, uneg, vpos, vneg){
  # This gets a single factor. Do MultiSMD to get multiple factors.
  # determines any missing values
  nas <- is.na(x)
  # accepts pre-calculated v matrix using CheckPMDV function
  v.init <- v
  # initialize x matrix to starting point
  xoo <- x
  # If there are any missing values, they are filled in with the mean value
  if(sum(nas)>0) xoo[nas] <- mean(x[!nas])
  # an initialized vector for v is generated from random x values
  oldv <- rnorm(ncol(x))
  # the number of iterations for convergence is specified as a parameter
  #   the default is set to be 20
  for(iter in 1:niter){
    # the criterion for convergence is a difference less than 1e-7
    if(sum(abs(oldv-v))>1e-7){
      # old v vector is replaced with updated v
      oldv <- v
      if(trace) cat(iter,fill=F)
      # update u #
      argu <- xoo%*%v
      if(upos) argu <- pmax(argu,0)
      if(uneg) argu <- pmin(argu,0)
      # Use method described in Sharifzadeh
      u <- matrix(argu/l2n(argu),ncol=1)
      # done updating u #
      # update v #
      argv <- t(u)%*%xoo
      if(vpos) argv <- pmax(argv,0)
      if(vneg) argv <- pmin(argv,0)
      # calculate appropriate sparse penalty
      lamv <- BinarySearch(argv, sumabsv)
      # perform soft threshold using sparseness penalty calculated using
      #   BinarySearch
      sv <- soft(argv,lamv)
      # convert v vector to matrix orm
      v <- matrix(sv/l2n(sv),ncol=1)
      # done updating v #
    }
  }
  # calculate d component
  d <- as.numeric(t(u)%*%(xoo%*%v))
  if(trace) cat(fill=TRUE)
  # return list of d, u, and v components along with initial v vector
  return(list(d=d, u=u, v=v, v.init=v.init))
}