#' A function that calculates the mean value of all non-missing values
#' of a vector
#' @param vec vector in which to calculate the mean of non-missing values
#' @keywords mean missing
#' @export
#' @examples
#' mean.na(vec)

# returns mean of a vector that consists of non-missing values
mean.na <- function(vec){
  return(mean(vec[!is.na(vec)]))
}