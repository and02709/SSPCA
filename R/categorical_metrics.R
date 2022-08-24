#' This function calculates the performance metrics for each categorical response.
#' @param ytrain training response dataset
#' @param ytest testing response dataset
#' @param y_hat fitted/predicted values
#' @param nresp number of response vectors present in datasets
#' @param resp.names gives list of names for each response
#' @keywords categorical metrics
#' @export
#' @examples categorical_metrics(ytrain, ytest, y_hat, nresp, resp.names)

categorical_metrics <- function(ytrain, ytest, y_hat, nresp, resp.names){
  y_hat_train <- y_hat[[1]]
  y_hat_test <- y_hat[[2]]
  
  y_train_class <- matrix(0, nrow=dim(as.matrix(ytrain))[1], ncol=dim(as.matrix(ytrain))[2])
  y_test_class <-  matrix(0, nrow=dim(as.matrix(ytest))[1], ncol=dim(as.matrix(ytest))[2])
  
  
  for(i in 1:nresp){
    if(nresp==1){
      y <- as.matrix(ytrain)
      y_hat <- as.matrix(y_hat_train)
      y_train_class <- (y==y_hat)
      acc_train <- mean(y_train_class)
      miscl_train <- 1-acc_train
    }
    else{
      y <- as.matrix(ytrain[,i])
      y_hat <- as.matrix(y_hat_train[,i])
      y_train_class[,i] <- (y==y_hat)
      acc_train <- colMeans(y_train_class)
      miscl_train <- 1-acc_train
    }
    
    
    if(nresp==1){
      y <- as.matrix(ytest)
      y_hat <- as.matrix(y_hat_test)
      y_test_class <- (y==y_hat)
      acc_test <- mean(y_test_class)
      miscl_test <- 1-acc_test
      
    }
    else{
      y <- as.matrix(ytest[,i])
      y_hat <- as.matrix(y_hat_test[,i])
      y_test_class[,i] <- (y==y_hat)
      acc_test <- colMeans(y_test_class)
      miscl_test <- 1-acc_test
    }
    
  }
  
  y_train_class <- as.matrix(y_train_class)
  y_test_class <- as.matrix(y_test_class)
  colnames(y_train_class) <- paste0("predicted ", resp.names)
  colnames(y_test_class) <- paste0("predicted ", resp.names)
  
  return(list(y_train_class=y_train_class, y_test_class=y_test_class, accuracy_train=acc_train, misclassification_error_train=miscl_train, accuracy_test=acc_test, misclassification_error_test=miscl_test))
  
  
}