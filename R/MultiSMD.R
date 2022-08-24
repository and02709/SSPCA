#' This function decomposes a matrix using singular value decomposition ad uses
#'   the SMD function.  This results in vectors corresponding to the loadings
#'   which have undergone shrinkage
#' @param x this is the matrix which is decomposed and the decomposition vectors
#'   corresponding to the loadings are shrunk
#' @param sumabsu this parameter determines by how much the decomposition vectors
#'   corresponding to the observations
#' @param sumabsv this parameter determines by how much the decomposition vectors
#'   corresponding to the loadings
#' @param K determines how many vectors should be calculated
#' @param niter number of iterations to allow for convergence
#' @param v allows a particular v matrix to be passed to the SMD function
#' @param trace option to deisplay potential warnings
#' @param upos allows the option to force the vector corresponding to 
#'   observations to correspond to the parallel maxima using pmax
#' @param uneg allows the option to force the vector corresponding to 
#'   observations to correspond to the parallel minima using pmin
#' @param vpos allows the option to force the vector corresponding to 
#'   loadings to correspond to the parallel maxima using pmax
#' @param vneg allows the option to force the vector corresponding to 
#'   loadings to correspond to the parallel minima using pmin
#' @keywords multiple singular matrix
#' @export
#' @examples MultiSMD(x, sumabsu, sumabsv, K=3, niter=20,v, trace=TRUE, upos, uneg, vpos, vneg)


# Function accepts the matrix x to decompose and employs sparseness using
#   L1 penalty on the vectors corresponding the loadings
MultiSMD <- function(x, sumabsu, sumabsv, K=3, niter=20,v, trace=TRUE, upos, uneg, vpos, vneg){
  # this identifies any missing values in the x matrix
  nas <- is.na(x)
  # initialize v to the matrix provided previously
  v.init <- v
  # pass x matrix to decompose
  xuse <- x
  # identify the number of singular values desired
  ds <- numeric(K)
  # initialize u and v matrices using normal values
  us <- matrix(0,nrow=nrow(x),ncol=K)
  vs <- matrix(0,nrow=ncol(x),ncol=K)
  # decompose x matrix into desired number of vectors
  for(k in 1:K){
    # calculate u,v, and singular values using SMD
    out <- SMD(xuse, sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,v=matrix(v[,k],ncol=1), trace=trace, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg)
    us[,k] <- out$u
    vs[,k] <- out$v
    ds[k] <- out$d
    # substract corresponding calculated values from original to achieve orthogonality
    res <- xuse - out$d*out$u%*%t(out$v)
    xuse[!nas] <- res[!nas] # [!nas] is new on July 24 2009
  }
  # return u, v, and singular values along with initialize v matrix
  return(list(u=us,v=vs,d=ds, v.init=v.init))
}