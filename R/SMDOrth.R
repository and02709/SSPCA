#' This function decomposes a matrix using singular value decomposition but 
#'   constructs the new u vectors such that they are orthogonal to the previous
#'   vectors
#' The resulting eigenvector corresponding to the loadings in then shrunk
#'   using an L1 penalty where the u vectors are orthogonal to previous u
#'   vectors.
#' @param x this is the matrix which is decomposed and the decomposition vectors
#'   corresponding to the loadings are shrunk
#' @param us this tells the us the previous u vectors for which we wish
#'   to construct new u vectors that are orthogonal
#' @param sumabsv this parameter determines by how much the decomposition vectors
#'   corresponding to the loadings
#' @param niter number of iterations to allow for convergence
#' @param trace option to deisplay potential warnings
#' @param v allows a particular loadings vector to be used
#' @param vpos allows the option to force the vector corresponding to 
#'   loadings to correspond to the parallel maxima using pmax
#' @param vneg allows the option to force the vector corresponding to 
#'   loadings to correspond to the parallel minima using pmin
#' @keywords singular matrix orthogonal
#' @export
#' @examples SMDOrth(x, us, sumabsv=NULL, niter=20, trace=TRUE,v, vpos, vneg)

# Function accepts the matrix x to decompose and employs sparseness using
#   L1 penalty on the vectors corresponding the loadings
SMDOrth <- function(x, us, sumabsv=NULL, niter=20, trace=TRUE,v, vpos, vneg){
  # Gets the next u for sparse PCA, using Trevor's method of requiring this new u to be orthog
  # to all previous us (contained in us)
  # determines which values of x may be missing
  nas <- is.na(x)
  # intitializes v vector to the one provided
  v.init <- v
  # initializes x matrix
  xoo <- x
  # any missing values are filled in with the mean value
  if(sum(nas)>0) xoo[nas] <- mean(x[!nas])
  # an initial u vector was generated using random normal values
  u <- rnorm(nrow(x))
  # an old u vector is generated again.  The assumption is that there will
  #   be enough difference between u and oldu to begin the algorithm
  oldu <- rnorm(nrow(x))
  # the old v vector is initialized using random normals
  oldv <- rnorm(ncol(x))
  # niter is a parameter but the default setting is 20
  for(iter in 1:niter){
    # the process will iterate if differences in either u or v are greater than
    # the value 1e-6
    if(sum(abs(oldu-u))>1e-6 || sum(abs(oldv-v))>1e-6){
      # update u and v
      oldu <- u
      oldv <- v
      if(trace) cat(iter,fill=F)
      # update u #
      argu <- xoo%*%v
      # use the lsfit to find the space orthogonal to the previous u vectors
      numer <- lsfit(y=argu, x=us, intercept=FALSE)$res
      # normalize the u vector so that it has norm 1
      u <- numer/l2n(numer)
      # done updating u #
      # update v #
      argv <- t(u)%*%xoo
      if(vpos) argv <- pmax(argv,0)
      if(vneg) argv <- pmin(argv,0)
      # caclulate sparsness shrinkage using BinarySearch
      lamv <- BinarySearch(argv, sumabsv)
      # perform soft thresholding on v vector
      sv <- soft(argv,lamv)
      # normalize v vector
      v <- matrix(sv/l2n(sv),ncol=1)
      # done updating v #
    }
  }
  # calculate d component
  d <- as.numeric(t(u)%*%(xoo%*%v))
  if(trace) cat(fill=TRUE)
  # return components d, u, v, and also initial v vector for calculations
  return(list(d=d,u=u,v=v, v.init=v.init))
}