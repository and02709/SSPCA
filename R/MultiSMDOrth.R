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
#' @keywords multiple singular matrix orthogonal
#' @export
#' @examples MultiSMDOrth(x,sumabsu,sumabsv,niter,K, trace, v, vpos, vneg)

# Function accepts the matrix x to decompose and employs sparseness using
#   L1 penalty on the vectors corresponding the loadings
MultiSMDOrth <- function(x,sumabsu,sumabsv,niter,K, trace, v, vpos, vneg){
  # checks sumabs for u to make sure shrinkage can occur on u 
  #   this has largely been disabled in this algorithm
  if(sumabsu < sqrt(nrow(x))){
    warning("sumabsu was less than sqrt(nrow(x)), so the orthogonal option was not implemented.")
    # This may create problems because I removed the U sparsity component
    return(MultiSMD(x,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,K=K, trace=trace, v=v, vpos=vpos, vneg=vneg))
  }
  # identfies potential missing values in the x matrix
  nas <- is.na(x)
  # accepts proposed v matrix corresponding to the loadings
  v.init <- v
  # assign x matrix to algorithm to use
  xuse <- x
  # generate initial singular values, u and v vectors
  ds <- numeric(K)
  us <- matrix(0,nrow=nrow(x),ncol=K)
  vs <- matrix(0,nrow=ncol(x),ncol=K)
  # begin the process of constructing sparse vecctors
  for(k in 1:K){
    # If no other eigenvectors are present, there is no need to orthogonalize
    #   the space to previous u matrices.  Therefore just use SMD
    if(k==1) out <- SMD(xuse,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,v=matrix(v[,k],ncol=1), trace=trace, upos=FALSE, uneg=FALSE,  vpos=vpos, vneg=vneg)
    # If there are other u vectors present, the space must be orthogonalized
    #   In this case use SMDOrth
    if(k>1) out <- SMDOrth(xuse,us[,1:(k-1)],sumabsv=sumabsv,niter=niter,v=matrix(v[,k],ncol=1), trace=trace,  vpos=vpos, vneg=vneg)
    us[,k] <- out$u
    vs[,k] <- out$v
    ds[k] <- out$d
    # I added this based on the MultiSMD function  It subtracts the components
    #   already identified to orthogonalize the space
    res <- xuse - out$d*out$u%*%t(out$v)
    # I added this  based on the MultiSMD function
    xuse[!nas] <- res[!nas] # [!nas] is new on July 24 2009
    
  }
  # returns the u matrix, v matrix, singular values, and initial v matrix
  return(list(u=us,v=vs,d=ds, v.init=v.init))
}