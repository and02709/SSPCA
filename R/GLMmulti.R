#' This function calculates predictions for each multinomial response using all
#'   eigenvectors together.
#' @param ytrain training response dataset
#' @param ztrain encoded training data
#' @param ytest testing response dataset
#' @param ztest encoded testing data
#' @param nresp number of response vectors present in datasets
#' @param resp.names gives list of names for each response
#' @keywords GLM multinomial
#' @export
#' @examples GLMmulti(ytrain, ztrain, ytest, ztest, nresp, resp.names)

GLMmulti <- function(ytrain, ztrain, ytest, ztest, nresp, resp.names){
  
  y_hat_train <- matrix(0, nrow=dim(as.matrix(ytrain))[1], ncol=dim(as.matrix(ytrain))[2])
  y_hat_test <- matrix(0, nrow=dim(as.matrix(ytest))[1], ncol=dim(as.matrix(ytest))[2])
  
  for(i in 1:nresp){
    if(nresp==1){
      y <- as.matrix(ytrain)
      npred <- dim(as.matrix(ztrain))[2]
      dftrain <- data.frame(y, ztrain)
      dftrain <- data.frame(y, ztrain)
      colnames(dftrain) <- c("y", paste0("z",1:npred))
      m1 <- multinom(y ~ ., data=dftrain)
      yhat <- predict(m1, dftrain)
      y_hat_train <- yhat
    }
    else{
      y <- as.matrix(ytrain[,i])
      npred <- dim(as.matrix(ztrain))[2]
      dftrain <- data.frame(y, ztrain)
      colnames(dftrain) <- c("y", paste0("z",1:npred))
      m1 <- multinom(y ~ ., data=dftrain)
      yhat <- predict(m1, dftrain)
      y_hat_train[,i] <- yhat
    }
    
    if(nresp==1){
      y <- as.matrix(ytest)
      dftest <- data.frame(y, ztest)
      colnames(dftest) <- c("y", paste0("z",1:npred))
      yhat <- predict(m1, dftest)
      y_hat_test <- yhat
    }
    else{
      y <- as.matrix(ytest[,i])
      dftest <- data.frame(y, ztest)
      colnames(dftest) <- c("y", paste0("z",1:npred))
      yhat <- predict(m1, dftest)
      y_hat_test[,i] <- yhat
    }
    
  }
  
  y_hat_train <- as.matrix(y_hat_train)
  y_hat_test <- as.matrix(y_hat_test)
  
  colnames(y_hat_train) <- paste0("hat ", resp.names)
  colnames(y_hat_test) <- paste0("hat ", resp.names)
  
  return(list(yhat_train=y_hat_train, yhat_test=y_hat_test))
  
}