#' This extracts the desired number of eigenvectors from the X matrix
#' Performs svd on x matrix in one of two ways.  First is if ncol > nrows,
#' svd performed on x t(x), otherwise svd performed on t(x) x.  Also, if there
#' are too many missing values, svd is performed on just x
#' @param v proposed eigenvectors.  Usually set to NULL
#' @param x matrix to extract eigenvectors
#' @param K number of desired eigenvectors
#' @keywords penalized matrix decomposition
#' @export
#' @examples
#' CheckPMDV(v,x,K)


# Normally accepts the v eigenvector of loadings, x input matrix, and number
# of desired eigenvalues.  In the case of SSPC, v is always null, and
# x represents the Psi matrix used to decompose the HSIC criterion.  K still
# represents the number of desired eigenvectors.
CheckPMDV <-  function(v,x,K){
  # this checks whether the vector v is non-null, a matrix, and contains more
  # than or equal to the number of desired columns
  if(!is.null(v) && is.matrix(v) && ncol(v)>=K){
    # select the desird number of columns
    v <- matrix(v[,1:K], ncol=K)
    
    # if the number of columns is greater than the number of rows in the desired
    # decomposition space
  } else if(ncol(x)>nrow(x)){
    # fill in missing values with means of remaining values
    x[is.na(x)] <- mean.na(x)
    # perform svd on deisred x by t(x) matrix and select K columns
    v <- matrix(t(x)%*%(safesvd(x%*%t(x))$v[,1:K]),ncol=K)
    # if there are missing values in v, perform svd on just x and extract
    #   K eigenvectors
    if(sum(is.na(v))>0) v <- matrix(safesvd(x)$v[,1:K], ncol=K)
    # This standardizes the v matrix so that all columns have length 1
    v <- sweep(v,2,apply(v, 2, l2n), "/")
    # If there are remaining missing values, halt algorithm
    if(sum(is.na(v))>0) stop("some are NA")
    # alternative method to extract eigenvectors if there are fewer cols than rows
  } else if (ncol(x)<=nrow(x)){
    x[is.na(x)] <- mean.na(x)
    v <- matrix(safesvd(t(x)%*%x)$v[,1:K],ncol=K)
  }
  return(v)
}