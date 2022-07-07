GLMbinSep <- function(ytrain, ytest, Z, nresp, resp.names){
  
  y_hat_train <- matrix(0, nrow=dim(as.matrix(ytrain))[1], ncol=dim(as.matrix(ytrain))[2])
  y_hat_test <- matrix(0, nrow=dim(as.matrix(ytest))[1], ncol=dim(as.matrix(ytest))[2])
  
  for(i in 1:nresp){
    if(nresp==1){
      y <- as.matrix(ytrain)
      ztrain <- Z[[i]]$Ztrain
      npred <- dim(as.matrix(ztrain))[2]
      dftrain <- data.frame(y, ztrain)
      colnames(dftrain) <- c("y", paste0("z",1:npred))
      m1 <- glm(y ~ .,family=binomial(link=logit), data=dftrain)
      yhatnum <- predict(m1, dftrain, type= "response")
      yhatchar <- as.character(ifelse(yhatnum < 0.5,0,1))
      yhatclass <- factor(yhatchar,levels=c("0","1"))
      y_hat_train <- yhatclass
    }
    else{
      y <- as.matrix(ytrain[,i])
      ztrain <- Z[[i]]$Ztrain
      npred <- dim(as.matrix(ztrain))[2]
      dftrain <- data.frame(y, ztrain)
      colnames(dftrain) <- c("y", paste0("z",1:npred))
      m1 <- glm(y ~ .,family=binomial(link=logit), data=dftrain)
      yhatnum <- predict(m1, dftrain, type= "response")
      yhatchar <- as.character(ifelse(yhatnum < 0.5,0,1))
      yhatclass <- factor(yhatchar,levels=c("0","1"))
      y_hat_train[,i] <- yhatclass
    }
    
    if(nresp==1){
      y <- as.matrix(ytest)
      ztrain <- Z[[i]]$Ztest
      npred <- dim(as.matrix(ztest))[2]
      dftest <- data.frame(y, ztest)
      colnames(dftest) <- c("y", paste0("z",1:npred))
      yhatnum <- predict(m1, dftest, type= "response")
      yhatchar <- as.character(ifelse(yhatnum < 0.5,0,1))
      yhatclass <- factor(yhatchar,levels=c("0","1"))
      y_hat_test[,i] <- yhatclass
    }
    else{
      y <- as.matrix(ytest[,i])
      ztrain <- Z[[i]]$Ztest
      npred <- dim(as.matrix(ztest))[2]
      dftest <- data.frame(y, ztest)
      colnames(dftest) <- c("y", paste0("z",1:npred))
      yhatnum <- predict(m1, dftest, type= "response")
      yhatchar <- as.character(ifelse(yhatnum < 0.5,0,1))
      yhatclass <- factor(yhatchar,levels=c("0","1"))
      y_hat_test <- yhatclass
    }
    
  }
  
  y_hat_train <- as.matrix(y_hat_train)
  y_hat_test <- as.matrix(y_hat_test)
  
  colnames(y_hat_train) <- paste0("hat ", resp.names)
  colnames(y_hat_test) <- paste0("hat ", resp.names)
  
  return(list(yhat_train=y_hat_train, yhat_test=y_hat_test))
}