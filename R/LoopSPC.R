#' This function cycles through each response and computes supervised principal
#'     component analysis eigenvectors.  This is caused by doing separate analysis.
#' @param ytrain training response dataset
#' @param xtrain training predictor dataset
#' @param ytest testing response dataset
#' @param xtest testing predictor dataset
#' @param ycont flag indicates response data is continuous
#' @param nresp number of response vectors present in datasets
#' @param strictEV flag whether to impose strict requiresments to accept number
#'   of eigenevectors or to select a more appropriate number
#' @param nTopEvecs number of desired eigenevectors
#' @keywords individual response supervised principal component analysis
#' @export
#' @examples LoopSPC(ytrain, xtrain, ytest, xtest, ycont, nresp, strictEV, nTopEvecs)

LoopSPC <- function(ytrain, xtrain, ytest, xtest, ycont, nresp, strictEV, nTopEvecs){
  
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
    
    
    
    # generate identity matrix that has dimensions equal to the training set
    I <- diag(nrow=nobs)
    # generate vector of 1's with length of the training set
    e <- rep(1, nobs)
    # generate square matrix of 1's with dimension of the training set
    H <- I -(1/nobs)*e%*%t(e)
    # generate the X matrix
    X <- as.matrix(xtrain)
    # generate centered matrix of kernel structure coupled with response matrix
    Q <- t(X)%*%H%*%L%*%H%*%X
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
    
    if(j==1){
      # Encode Training Data
      Ztrain <- t(t(Uvectors)%*%t(X))
      # Encode Testing Data
      Ztest <- t(t(Uvectors)%*%t(as.matrix(xtest)))
      temp_hold <- list(Ztrain=Ztrain, Ztest=Ztest, U=Uvectors, lambda=Uval)
      #assign(name_list[j], temp_hold)
      list_of_list <- list(temp_hold)
    }
    else{
      ZtrTemp <- t(t(Uvectors)%*%t(X))
      ZtsTemp <- t(t(Uvectors)%*%t(as.matrix(xtest)))
      
      # Ztrain <- cbind(Ztrain, ZtrTemp)
      # Ztest <- cbind(Ztest, ZtsTemp)
      temp_hold <- list(Ztrain=ZtrTemp, Ztest=ZtsTemp, U=Uvectors, lambda=Uval)
      #assign(name_list[j], temp_hold)
      list_of_list <- append(list_of_list, list(temp_hold))
    }
    
    
  }
  
  return(list_of_list)
  
}