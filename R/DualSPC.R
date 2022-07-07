

DualSPC <- function(ytrain, xtrain, ytest, xtest, ycont, ybinary, nresp, strictEV, nTopEvecs){
  
  nobs <- nrow(as.matrix(ytrain))
  
  if(ycont ){
    L <- as.matrix(ytrain)%*%t(as.matrix(ytrain))
  }
  
  if(!ycont && ybinary){
    
    yf <- recPack(ytrain=ytrain)
    yf <- as.factor(yf)
    L<-matrix(0,nobs,nobs)
    for (i in levels(yf)){
      tmp<-yf==i
      L[tmp,tmp]<-1
    }
  }
  
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
  
  
  # Eigen decomposition of matrix of kernel response matrix L
  Eigendecomp <- eigen(L)
  # Need to retain eigenvectors to reconstruct Delta
  U <- Eigendecomp$vectors
  # Eigenvalues of kernel response matrix L
  EV <- Eigendecomp$values
  # Generation of diagonal matrix of square root eigenvalues of L
  Sigmat <- diag(sqrt(zapsmall(EV)))
  # Generation of Delta such that t(Delta)%*%Delta=L
  Delta <- zapsmall(U%*%Sigmat%*%t(U))
  # generate identity matrix that has dimensions equal to the training set
  I <- diag(nrow=nobs)
  # generate vector of 1's with length of the training set
  e <- rep(1, nobs)
  # generate square matrix of 1's with dimension of the training set
  H <- I -(1/nobs)*e%*%t(e)
  # generate the X matrix
  X <- as.matrix(xtrain)
  # Generation of Psi matrix
  Psi <- t(Delta)%*%H%*%X
  # generate Dual Supervised form
  Q <- Psi%*%t(Psi)
  # Eigen decomposition of Q matrix
  Utotal <- eigen(Q, only.values=F)
  # This object stores the eigenvalues
  Uval <- Utotal$values
  # This object stores the eigenvectors
  Uvec <- Utotal$vectors
  
  # Some of the eigenvalues and eigenvectors may be complex if the
  #   resulting Q matrix has values approaching zero.  Therefore, 
  #   we must remove the complex component or future operations will not work
  if(is.complex(Uval)){
    Uval <- Re(Uval)
    Uvec <- Re(Uvec)
  }
  
  if(strictEV){
    #zapUval <- zapsmall(Uval)
    Uvalues <- Uval[1:nTopEvecs]
    Uvectors <- Uvec[,1:nTopEvecs]
    numvec <- length(Uvalues)
    
  }
  
  else{
    zapUval <- zapsmall(Uval)
    ntake <- sum(zapUval > 0)
    numvec <- 0
    
    if(ntake==1 || ntake < nTopEvecs){
      Indec <- which(zapUval>0)
      numvec <- length(Indec)
      
      if(numvec>0){
        Uvalues <- zapUval[Indec]
        Uvectors <- Uvec[,Indec]
      }
    }
    
    else{
      Uvalues <- zapUval[1:nTopEvecs]
      Uvectors <- Uvec[,1:nTopEvecs]
      numvec <- length(Uvalues)
    }
    
  }
  
  
  # This code addresses the circumstance that no eigenvectors are found
  if(numvec==0) stop("No eigenvalues greater than 0")
  
  # Generate Encoded Training Data
  if(numvec==1){
    Siggy <- sqrt(Uvalues)
    Ut <- t(t(Psi)%*%Uvectors%*%t(solve(Siggy)))
  }
  else{
    Siggy <- diag(sqrt(Uvalues))
    Ut <- t(t(Psi)%*%Uvectors%*%t(solve(Siggy)))
  }
  Ztrain <- t(Ut%*%t(X))
  # Generate Encoded Testing Data
  Ztest <- t(Ut%*%t(as.matrix(xtest)))
  
  
  return(list(Z=Ztrain, z=Ztest, U=t(Ut), lambda=Uval))
}
