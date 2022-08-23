# This function assigns a category based on binary response permutation
# @param ytrain accepts matrix of response variables
# @keywords
# @export
# @examples
# recPack(ytrain)

# Function accepts ytrain matrix of response variables
recPack <- function(ytrain){
  # converts ytrain data to a matrix
  ytemp <- as.matrix(ytrain)
  # number of observations
  n <- nrow(ytemp)
  # number of response variables
  p <- ncol(ytemp)
  # builds matrix of possible response permutations
  Key <- t(as.matrix(PermBinary(p)))
  # the following loop converts each response into a numerical value of 1 and 0
  #   depending on the number of responses the number of integers is equal to p
  #   the conversion to base ten allows the uniqueness to remain
  for(i in 1:p){
    if(i==1){
      yt <- (10^(p-i))*ytemp[,i]
    }
    else{
      yt <- yt + (10^(p-i))*ytemp[,i]
    }
  }
  # This loop performs the same conversion on the key
  for(i in 1:p){
    if(i==1){
      kt <- (10^(p-i))*Key[,i]
    }
    else{
      kt <- kt + (10^(p-i))*Key[,i]
    }
  }
  
  j <- nrow(Key)
  y <- rep(0, n)
  # This loop assigns categories based on the numerical values generated in the
  #   previous two loops
  for(i in 1:j){
    y[which(yt==kt[i])] <- i
  }
  return(y)
}