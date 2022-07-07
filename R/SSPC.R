# This function supervises the decomposition of the matrix X into vectors
#   corresponding to the observations and sparse vectors corresponding to 
#   the loadings.  It handles all operations with getting the desired
#   eigenvectors, and encodes the data.  This information is passed
#   to subsequent functions
# @param xtrain training predictor dataset
# @param ytrain training response dataset
# @param xtest testing predictor dataset
# @param ycont flag for whether the response data is continuous
# @param ybinary flag for whether the response data is binary
# @param nresp number of response vectors in the training response dataset
# @param sumabsv this parameter determines by how much the decomposition vectors
#   corresponding to the loadings
# @param niter number of iterations to allow for convergence
# @param K determines how many vectors should be calculated
# @param v allows a particular v matrix to be passed to PMDL1
# @param trace option to deisplay potential warnings
# @param center this tells the algorithm whether centering is required
#   this has been handled in other programs so has been disabled
# @param cnames gives potential column names
# @param upos allows the option to force the vector corresponding to 
#   observations to correspond to the parallel maxima using pmax
# @param uneg allows the option to force the vector corresponding to 
#   observations to correspond to the parallel minima using pmin
# @param vpos allows the option to force the vector corresponding to 
#   loadings to correspond to the parallel maxima using pmax
# @param vneg allows the option to force the vector corresponding to 
#   loadings to correspond to the parallel minima using pmin
# @param compute.pve no longer useful
# @param strictEV determines whether the specified K must be strictly followed
#   or whether a more appropriate number of vectors can be used

# Function accepts the training and testing predictor datasets and the
#   training response datasets
SSPC <- function(xtrain, ytrain, xtest, ycont, ybinary, nresp, sumabsv=4, niter=20, K=1, orth=TRUE, trace=TRUE, v=NULL, center=FALSE, cnames=NULL, vpos=FALSE, vneg=FALSE, compute.pve=TRUE, strictEV=TRUE){
  if(vpos&&vneg) stop("Cannot constrain elements to be positive AND negative.")
  # calculates the nunber of observations present in the training dataset
  nobs <- nrow(as.matrix(ytrain))
  
  # if response is continuous, use the yy^t kernel
  if(ycont){
    L <- as.matrix(ytrain)%*%t(as.matrix(ytrain))
  }
  
  # if response is binary, use the delta kernel
  if(!ycont && ybinary){
    
    yf <- recPack(ytrain=ytrain)
    yf <- as.factor(yf)
    L<-matrix(0,nobs,nobs)
    for (i in levels(yf)){
      tmp<-yf==i
      L[tmp,tmp]<-1
    }
  }
  
  # this produces a kernel matrix for categorical data, but only one vector
  #   cannot handle more than one categories
  if(!ycont && !ybinary){
    if(nresp>1) stop("Cannot transform multivariable responses with categorical variables into response kernel")
    yf <- as.factor(ytrain)
    L<-matrix(0,nobs,nobs)
    for (i in levels(yf)){
      tmp<-yf==i
      L[tmp,tmp]<-1
    }
  }
  
  #if(improper.Y.Binary && !ycont && ybinary){
  #  L <- as.matrix(ytrain)%*%t(as.matrix(ytrain))
  #}
  
  # Decompose the response kernel and reconstruct it such that you
  #   get a delta matrix that is the 1/2 matrix of the L kernel
  Eigendecomp <- eigen(L)
  U <- Eigendecomp$vectors
  EV <- Eigendecomp$values
  Sigmat <- diag(sqrt(zapsmall(EV)))
  Delta <- zapsmall(U%*%Sigmat%*%t(U))
  
  # centering matrix
  H <- diag(1, nrow(as.matrix(ytrain))) - 1/nrow(as.matrix(ytrain))*rep(1, nrow(as.matrix(ytrain)))%*%t(rep(1, nrow(as.matrix(ytrain))))
  
  # HSIC matrix to be decomposed
  Psi <- t(Delta)%*%H%*%as.matrix(xtrain)
  
 
  
  Q <- Psi%*%t(Psi)
  EigenValues <- eigen(Q)[[1]]

  # Provided strictEV is set to FALSE, this allows a more reasonable number
  #   of eigenvectors to be selected
  if(!strictEV){
    EigenValues2 <- zapsmall(eigen(Q)[[1]])
    K <- sum(EigenValues2 > 0)
  }
  
  
  # calls PMDL1 to generate sparse solutions
  out <- PMDL1(Psi,sumabsu=sqrt(nrow(xtrain)), sumabsv=sumabsv, niter=niter,K=K,orth=orth,trace=trace,v=v,center=center,cnames=cnames, upos=FALSE, uneg=FALSE, vpos=vpos, vneg=vneg)
  
  # sparse eigenvector loadings
  V <- as.matrix(out$v)
  
  #encode training and testing data
  ztrain <- as.matrix(xtrain)%*%V
  ztest <- as.matrix(xtest)%*%V
  
  #class(out) <- "SPC"
  return(list(Z=ztrain, z=ztest, U=V, lambda=EigenValues))
}