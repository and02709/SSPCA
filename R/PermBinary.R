# This function caclulates all possible response permutations for a series
# of binary response variables
# @param n the number of binary responses
# @keywords
# @export
# @examples
# PermBinary(n)

# Function accepts n, the number of binary responses
PermBinary <- function(n){ 
  # accounts for the total number of permutations
  n.perms <- 2^n 
  # establishes array that contains each possible response permutation
  array <- matrix(0,nrow=n,ncol=n.perms) 
  # array <- big.matrix(n, n.perms, type='integer', init=-5) 
  for(i in 1:n){ 
    div.length <- ncol(array)/(2^i) 
    div.num <- ncol(array)/div.length 
    end <- 0 
    while(end!=ncol(array)){ 
      end <- end +1 
      start <- end + div.length 
      end <- start + div.length -1 
      array[i,start:end] <- 1 
    } 
  } 
  return(array) 
} 