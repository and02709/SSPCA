#' A soft threshold function
#' This function shrinks a vector by a given amount
#' @param x  This parameter accepts a vector to be shrunk
#' @param d  This parameter represents the magnitude fo the shrinkage on vector x
#' @keywords soft threshold
#' @export
#' @examples
#' soft(x,d)


soft <- function(x,d){
  return(sign(x)*pmax(0, abs(x)-d))
}