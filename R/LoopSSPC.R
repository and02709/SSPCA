#' This function cycles through each response and computes supervised principal
#'     component analysis eigenvectors.  This is caused by doing separate analysis.
#'     This invokes the SSPCA function to compute sparse eigenevectors
#'     for each response.
#' @param xtrain training predictor dataset
#' @param ytrain training response dataset
#' @param xtest testing predictor dataset
#' @param ycont flag indicates response data is continuous
#' @param nresp number of response vectors present in datasets
#' @param sumabsv causes shrinkage to occur to the loadings
#' @param niter the number of iterations allowed to achieve convergence
#' @param K number of eigenevectors allowed to be calculated
#' @param orth determines whether orthogonal vectors are required
#' @param trace this allows debugging
#' @param v focuses on the loadings vector
#' @param center this centers the predictor columns
#' @param cnames this gives column names
#' @param vpos this determines whether the matrix is to be positive only
#' @param vneg this determines whether the matrix is to be negative only
#' @param compute.pve computes cv using missing values
#' @param strictEV flag whether to impose strict requiresments to accept number
#'   of eigenevectors or to select a more appropriate number
#' @keywords indiviual sparse supervised principal component analysis
#' @export
#' @examples LoopSSPC(xtrain, ytrain, xtest, ycont, nresp, sumabsv=4, niter=20, 
#'   K=1, orth=TRUE, trace=TRUE, v=NULL, center=FALSE, 
#'   cnames=NULL, vpos=FALSE, vneg=FALSE, compute.pve=TRUE, 
#'   strictEV=TRUE)

LoopSSPC <- function(xtrain, ytrain, xtest, ycont, nresp, sumabsv=4, niter=20, 
                     K=1, orth=TRUE, trace=TRUE, v=NULL, center=FALSE, 
                     cnames=NULL, vpos=FALSE, vneg=FALSE, compute.pve=TRUE, 
                     strictEV=TRUE){
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