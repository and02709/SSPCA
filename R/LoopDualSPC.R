

LoopDualSPC <- function(ytrain, xtrain, ytest, xtest, ycont, nresp, strictEV, nTopEvecs){
  
  nobs <- nrow(ytrain)
  
  #name_list <- paste0("list_", 1:nresp)
  
  for(j in 1:nresp){
    
    y <- ytrain[,j]
    
    if(ycont){
      L <- as.matrix(y)%*%t(as.matrix(y))
    }
    else{
      yf <- as.factor(y)
      L<-matrix(0,nobs,nobs)
      for (i in levels(yf)){
        tmp<-yf==i
        L[tmp,tmp]<-1
      }
    }
    
    
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
    I <- diag(nrow=nrow(ytrain))
    # generate vector of 1's with length of the training set
    e <- rep(1, nrow(ytrain))
    # generate square matrix of 1's with dimension of the training set
    H <- I -(1/nrow(ytrain))*e%*%t(e)
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
        
        if(ntake>0){
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
    
    if(numvec==1){
      Siggy <- sqrt(Uvalues)
      Ut <- t(t(Psi)%*%Uvectors%*%t(solve(Siggy)))
    }
    else{
      Siggy <- diag(sqrt(Uvalues))
      Ut <- t(t(Psi)%*%Uvectors%*%t(solve(Siggy)))
    }
    #Ztrain <- t(Ut%*%t(X))
    # Generate Encoded Testing Data
    #Ztest <- t(Ut%*%t(as.matrix(xtest)))
    
    if(j==1){
      # Encode Training Data
      Ztrain <- t(Ut%*%t(X))
      # Encode Testing Data
      Ztest <- t(Ut%*%t(as.matrix(xtest)))
      temp_hold <- list(Ztrain=Ztrain, Ztest=Ztest, U=t(Ut), lambda=Uval)
      #assign(name_list[j], temp_hold)
      list_of_list <- list(temp_hold)
    }
    else{
      ZtrTemp <- t(Ut%*%t(X))
      ZtsTemp <- t(Ut%*%t(as.matrix(xtest)))
      
      #Ztrain <- cbind(Ztrain, ZtrTemp)
      #Ztest <- cbind(Ztest, ZtsTemp)
      temp_hold <- list(Ztrain=ZtrTemp, Ztest=ZtsTemp, U=t(Ut), lambda=Uval)
      #assign(name_list[j], temp_hold)
      list_of_list <- append(list_of_list, list(temp_hold))
    }
    
    
  }
  
  
  return(list_of_list)
  
}
