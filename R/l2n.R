#' A function calculates the l2 norm
#' This function calculates the l2n norm of a given vector "vec"
#' If the l2 norm is zero, this function then returns a value of 0.05 
#' so as not to result in a division by zero
#' @param vec Vector from which to calculate l2 norm
#' @keywords l2 norm
#' @export
#' @examples
#' l2n(vec)

# accepts vec parameter to calculate l2 norm
l2n <- function(vec){
  # calculates l2 norm
  a <- sqrt(sum(vec^2))
  # if l2 norm is 0, return 0.05 as to not be dividing by 0
  if(a==0) a <- .05
  # return l2 norm
  return(a)
}
