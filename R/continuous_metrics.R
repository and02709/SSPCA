#' This function calculates the performance metrics for each continuous response.
#' @param ytrain training response dataset
#' @param ytest testing response dataset
#' @param y_hat fitted/predicted values
#' @param nresp number of response vectors present in datasets
#' @param resp.names gives list of names for each response
#' @keywords 
#' @export
#' @examples continuous_metrics(ytrain, ytest, y_hat, nresp, resp.names)

continuous_metrics <- function(ytrain, ytest, y_hat, nresp, resp.names){
  y_hat_train <- y_hat[[1]]
  y_hat_test <- y_hat[[2]]
  
  y_train_SqEr <- matrix(0, nrow=dim(as.matrix(ytrain))[1], ncol=dim(as.matrix(ytrain))[2])
  y_test_SqEr <-  matrix(0, nrow=dim(as.matrix(ytest))[1], ncol=dim(as.matrix(ytest))[2])
  
  
  for(i in 1:nresp){
    if(nresp==1){
      y <- as.matrix(ytrain)
      y_hat <- as.matrix(y_hat_train)
      y_train_SqEr <- (y-y_hat)^2
      MSE_train <- mean(y_train_SqEr)
      RMSE_train <- sqrt(MSE_train)
    }
    else{
      y <- as.matrix(ytrain[,i])
      y_hat <- as.matrix(y_hat_train[,i])
      y_train_SqEr[,i] <- (y-y_hat)^2
      MSE_train <- colMeans(y_train_SqEr)
      RMSE_train <- sqrt(MSE_train)
    }
    
    
    if(nresp==1){
      y <- as.matrix(ytest)
      y_hat <- as.matrix(y_hat_test)
      y_test_SqEr <- (y-y_hat)^2
      MSE_test <- mean(y_test_SqEr)
      RMSE_test <- sqrt(MSE_test)
    }
    else{
      y <- as.matrix(ytest[,i])
      y_hat <- as.matrix(y_hat_test[,i])
      y_test_SqEr[,i] <- (y-y_hat)^2
      MSE_test <- colMeans(y_test_SqEr)
      RMSE_test <- sqrt(MSE_test)
    }
    
  }
  
  y_train_SqEr <- as.matrix(y_train_SqEr)
  y_test_SqEr <- as.matrix(y_test_SqEr)
  
  colnames(y_train_SqEr) <- paste0("Square Error ", resp.names)
  colnames(y_test_SqEr) <- paste0("Square Error ", resp.names)
  
  return(list(y_train_SqEr=y_train_SqEr, y_test_SqEr=y_test_SqEr, MSE_train=MSE_train, RMSE_train=RMSE_train, MSE_test=MSE_test, RMSE_test=RMSE_test))
  
  
}