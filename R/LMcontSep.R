#' This function calculates yhat for each continuous response separately.
#' @param ytrain training response dataset
#' @param ytest testing response dataset
#' @param Z encoded pca responses
#' @param nresp number of response vectors present in datasets
#' @param resp.names gives list of names for each response
#' @keywords linear model continuous separate
#' @export
#' @examples LMcontSep(ytrain, ytest, Z, nresp, resp.names)

LMcontSep <- function(ytrain, ytest, Z, nresp, resp.names){
  
  y_hat_train <- matrix(0, nrow=dim(as.matrix(ytrain))[1], ncol=dim(as.matrix(ytrain))[2])
  y_hat_test <- matrix(0, nrow=dim(as.matrix(ytest))[1], ncol=dim(as.matrix(ytest))[2])
  
  for(i in 1:nresp){
    if(nresp==1){
      y <- as.matrix(ytrain)
      ztrain <- Z[[i]]$Ztrain
      npred <- dim(as.matrix(ztrain))[2]
      dftrain <- data.frame(y, ztrain)
      colnames(dftrain) <- c("y", paste0("z",1:npred))
      m1 <- lm(y~., data=dftrain)
      yhat <- predict(m1, dftrain)
      y_hat_train <- yhat
    }
    else{
      y <- as.matrix(ytrain[,i])
      ztrain <- Z[[i]]$Ztrain
      npred <- dim(as.matrix(ztrain))[2]
      dftrain <- data.frame(y, ztrain)
      colnames(dftrain) <- c("y", paste0("z",1:npred))
      m1 <- lm(y~., data=dftrain)
      yhat <- predict(m1, dftrain)
      y_hat_train[,i] <- yhat
    }
    
    if(nresp==1){
      y <- as.matrix(ytest)
      ztest <- Z[[i]]$Ztest
      npred <- dim(as.matrix(ztest))[2]
      dftest <- data.frame(y, ztest)
      colnames(dftest) <- c("y", paste0("z",1:npred))
      yhat <- predict(m1, dftest)
      y_hat_test <- yhat
    }
    else{
      y <- as.matrix(ytest[,i])
      ztest <- Z[[i]]$Ztest
      npred <- dim(as.matrix(ztest))[2]
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