

LoopSSPC <- function(xtrain, ytrain, xtest, ycont, nresp, sumabsv=4, niter=20, K=1, orth=TRUE, trace=TRUE, v=NULL, center=FALSE, cnames=NULL, vpos=FALSE, vneg=FALSE, compute.pve=TRUE, strictEV=TRUE){
  if(vpos&&vneg) stop("Cannot constrain elements to be positive AND negative.")
  
  #name_list <- paste0("list_", 1:nresp)
  
  for(j in 1:nresp){
    
    y <- ytrain[,j]
    
    if(ycont){
      L <- y%*%t(y)
    }
    else{
      yf <- as.factor(y)
      L<-matrix(0,nrow(y),nrow(y))
      for (i in levels(yf)){
        tmp<-yf==i
        L[tmp,tmp]<-1
      }
    }
    
    Eigendecomp <- eigen(L)
    U <- Eigendecomp$vectors
    EV <- Eigendecomp$values
    Sigmat <- diag(sqrt(zapsmall(EV)))
    Delta <- zapsmall(U%*%Sigmat%*%t(U))
    
    H <- diag(1, nrow(as.matrix(y))) - 1/nrow(as.matrix(y))*rep(1, nrow(as.matrix(y)))%*%t(rep(1, nrow(as.matrix(y))))
    
    Psi <- t(Delta)%*%H%*%as.matrix(xtrain)
    
    Q <- Psi%*%t(Psi)
    EigenValues <- eigen(Q)[[1]]
    
    if(!strictEV){
      EigenValues2 <- zapsmall(eigen(Q)[[1]])
      K <- sum(EigenValues2 > 0)
    }
    
    
    out <- PMDL1(Psi,sumabsu=sqrt(nrow(xtrain)), sumabsv=sumabsv, niter=niter,K=K,orth=orth,trace=trace,v=v,center=center,cnames=cnames, upos=FALSE, uneg=FALSE, vpos=vpos, vneg=vneg)
    
    V <- as.matrix(out$v)
    
    if(j==1){
      # Encode Training Data
      Ztrain <- as.matrix(xtrain)%*%V
      # Encode Testing Data
      Ztest <- as.matrix(xtest)%*%V
      temp_hold <- list(Ztrain=Ztrain, Ztest=Ztest, U=V, lambda=EigenValues)
      #assign(name_list[j], temp_hold)
      list_of_list <- list(temp_hold)
    }
    else{
      ZtrTemp <- as.matrix(xtrain)%*%V
      ZtsTemp <- as.matrix(xtest)%*%V
      
      #Ztrain <- cbind(Ztrain, ZtrTemp)
      #Ztest <- cbind(Ztest, ZtsTemp)
      temp_hold <- list(Ztrain=ZtrTemp, Ztest=ZtsTemp, U=V, lambda=EigenValues)
      #assign(name_list[j], temp_hold)
      list_of_list <- append(list_of_list, list(temp_hold))
    }
    
    
  }
  
  
  return(list_of_list)
  
}