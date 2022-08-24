#' This function converts the eigenecvectors into a new form.
#' @param Z contains eigenvectors
#' @param resp.names gives list of names for each response
#' @param nresp number of response vectors present in datasets
#' @keywords 
#' @export
#' @examples restructuring.func(Z, resp.names, nresp)

restructuring.func <- function(Z, resp.names, nresp){
  if(length(Z)!=nresp) stop("Incorrect number of responses")
  
  
  nobs <- length(Z[[1]]$Ztrain)
  Ztrain <- matrix(0, nobs, nresp)
  for(i in 1:nresp){
    Ztrain[,i] <- Z[[i]]$Ztrain
  }
  colnames(Ztrain) <- paste0("Z ", resp.names)
  
  nobs <- length(Z[[1]]$Ztest)
  Ztest <- matrix(0, nobs, nresp)
  for(i in 1:nresp){
    Ztest[,i] <- Z[[i]]$Ztest
  }
  colnames(Ztest) <- paste0("z ", resp.names)
  
  nload <- length(Z[[1]]$U)
  U <- matrix(0, nload, nresp)
  for(i in 1:nresp){
    U[,i] <- Z[[i]]$U
  }
  colnames(U) <- paste0("U ", resp.names)
  
  nload <- length(Z[[1]]$lambda)
  Lambda <- matrix(0, nload, nresp)
  for(i in 1:nresp){
    Lambda[,i] <- Z[[i]]$lambda
  }
  colnames(Lambda) <- paste0("lambda ", resp.names)
  
  return(list(Z=Ztrain, z=Ztest, U=U, lambda=Lambda))
  
  
}