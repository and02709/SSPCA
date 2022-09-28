#' A function that calculates the maximum shrinkage for a given sparse penalty
#' This function examines te maximum elements of the eigenvector
#' and determines the level of shrinkage possible to meet the sparse penalty
#' @param argu This is the eigenvector to undergo l1 penalization
#' @param sumabs This is the sparse penalty 
#' @keywords shrinkage
#' @export
#' @examples
#' BinarySearch(argu,sumabs)

# The input parameters argu and sumabs that represent the eigenvector and
#   sparse penalty
BinarySearch <- function(argu,sumabs){
  # if the l2 norm is 0 or the shrinkage penalty is larger than the l1 norm
  #   of the vector no shrinkage will occur
  if(l2n(argu)==0 || sum(abs(argu/l2n(argu)))<=sumabs) return(0)
  # the limits for shrinkage are taken to be the minimum 0 and maximum 
  #   approximately the absolute value of the largest loading
  lam1 <- 0
  lam2 <- max(abs(argu))-1e-5
  iter <- 1
  # iteratively check whether the median sparse setting is greater or smaller
  #   than the l1 of the shrunk vector of interest denoted as su
  while(iter < 150){
    # shrink eigenvector by median of sparse settings
    su <- soft(argu,(lam1+lam2)/2)
    # if l1 norm is less than sparse penalty, increase the sparse settings
    if(sum(abs(su/l2n(su)))<sumabs){
      lam2 <- (lam1+lam2)/2
    } else {
    # if l1 norm is greater than sparse penalty, decrease sparse settings
      lam1 <- (lam1+lam2)/2
    }
    # return sparse settings if upper and lower sparse settings converge
    if((lam2-lam1)<1e-6) return((lam1+lam2)/2)
    iter <- iter+1
  }
  # warning if sparse penalties did not converge
  warning("Didn't quite converge")
  return((lam1+lam2)/2)
}

