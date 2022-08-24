#' This function supervises the decomposition of the matrix X into vectors
#'   corresponding to the observations and sparse vectors corresponding to 
#'   the loadings
#' @param x this is the matrix which is decomposed and the decomposition vectors
#'   corresponding to the loadings are shrunk
#' @param sumabsu this parameter determines by how much the decomposition vectors
#'   corresponding to the observations
#' @param sumabsv this parameter determines by how much the decomposition vectors
#'   corresponding to the loadings
#' @param niter number of iterations to allow for convergence
#' @param K determines how many vectors should be calculated
#' @param v allows a particular v matrix to be passed to the MultiSMD and
#'   MultiSMDOrth functions.  Its default is set to null because
#'   CheckPMDV is called in this function
#' @param trace option to deisplay potential warnings
#' @param orth specifies that MultiSMDOrth is to be called for any number of 
#'   vectors that is greater than 1
#' @param center this tells the algorithm whether centering is required
#'   this has been handled in other programs so has been disabled
#' @param rnames gives potential row names
#' @param cnames gives potential column names
#' @param upos allows the option to force the vector corresponding to 
#'   observations to correspond to the parallel maxima using pmax
#' @param uneg allows the option to force the vector corresponding to 
#'   observations to correspond to the parallel minima using pmin
#' @param vpos allows the option to force the vector corresponding to 
#'   loadings to correspond to the parallel maxima using pmax
#' @param vneg allows the option to force the vector corresponding to 
#'   loadings to correspond to the parallel minima using pmin
#' @keywords partial matrix decomposition
#' @export
#' @examples PMDL1(x,sumabs=.4,sumabsu=NULL,sumabsv=NULL,niter=20,
#'   K=1,v=NULL, trace=TRUE, orth=TRUE, center=TRUE, rnames=NULL, cnames=NULL, 
#'   upos=upos, uneg=uneg, vpos=vpos, vneg=vneg)


# Function accepts the matrix x to decompose and employs sparseness using
#   L1 penalty on the vectors corresponding the loadings
PMDL1 <- function(x,sumabs=.4,sumabsu=NULL,sumabsv=NULL,niter=20,K=1,v=NULL, trace=TRUE, orth=TRUE, center=TRUE, rnames=NULL, cnames=NULL, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg){
  
  # gives warning if sparseness is not possible due to numbers of elements
  if(orth){
    if(is.null(sumabsu) || sumabsu < sqrt(nrow(x))){
      orth <- FALSE
      warning("Orth option ignored because sparse PCA results only when sumabsu equals sqrt(nrow(x))")
    }
  }
  #  if((is.null(sumabsu) && !is.null(sumabsv)) || (is.null(sumabsv) && !is.null(sumabsu))) warning("Sumabsu and sumabsv BOTH must be input when type=standard. Since only one was given, it was ignored and sumabs was used instead.")
  # if no sumabs arguments given, it will fill this in based on the number
  #   of elements in rows or columns
  if(is.null(sumabsu) || is.null(sumabsv)){
    sumabsu <- sqrt(nrow(x))*sumabs
    sumabsv <- sqrt(ncol(x))*sumabs
  }
  call <-  match.call()
  # warning about centering, largely unneccessary at this point
  if(trace && (abs(mean.na(x)) > 1e-15)) warning("PMDL1 was run without first subtracting out the mean of x because this function does not perform the centering option by column.")
  # algoirthm is halted if sumabs values are too low or too high for the given
  #   dimensions of the matrix
  if(!is.null(sumabsu) && (sumabsu<1 || sumabsu>sqrt(nrow(x)))) stop("sumabsu must be between 1 and sqrt(n)")
  if(!is.null(sumabsv) && (sumabsv<1 || sumabsv>sqrt(ncol(x)))) stop("sumabsv must be between 1 and sqrt(p)")
  # This generates a reasonable initial starting point for decomposition
  #   steps
  v <- CheckPMDV(v,x,K)
  # This step is called only if no explicit instructions to ignore orthogonality
  #   are used
  if(K>1 && !orth) out <- (MultiSMD(x,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,K=K, trace=trace, v=v, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg))
  # This is the default option that is called that enures orthogonality 
  if(K>1 && orth) out <- MultiSMDOrth(x,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,K=K, trace=trace, v=v,  vpos=vpos, vneg=vneg)
  # called only if 1 vector is desired.  No need to orthogonlize
  if(K==1) out <- SMD(x,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter, trace=trace, v=v, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg)
  # returns all info in new object
  obj <- (list(u=out$u,v=out$v, d=out$d, v.init=out$v.init, call=call, sumabsu=sumabsu, sumabsv=sumabsv, rnames=rnames, cnames=cnames, K=K, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg))
  class(obj) <- "PMDL1"
  return(obj)
}