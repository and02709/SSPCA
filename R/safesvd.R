#' A function that performs SVD
#' This function iterates 10 times to attempt to decompose a matrix
#' through the singular value decomposition.  If the function is unable
#' to decompose the matrix after 10 tries, it will construct a random matrix
#' using rnorm with the same dimensions as the original matrix.  This random
#' matrix will be decomposed using the SVD and returned to the user.
#' @param x  This parameter accepts the matrix to undergo SVD
#' @keywords singular value decomposition
#' @export
#' @examples
#' safesvd(x)

safesvd <- function(x){
  # sets the index to iteratively attempt SVD
  i <- 1
  # attempt svd decomposition using try functiopn
  # this allows the function to recover if svd(x) fails
  out <- try(svd(x), silent=TRUE)
  # iteratively attempt svd using try function
  while(i<10 && class(out)=="try-error"){
    out <- try(svd(x), silent=TRUE)
    i <- i+1
  }
  # if svd on x fails, construct new matrix with same dimensions as x
  # conduct svd on new matrix
  if(class(out)=="try-error") out <- svd(matrix(rnorm(nrow(x)*ncol(x)), ncol=ncol(x)))
  return(out)
}
