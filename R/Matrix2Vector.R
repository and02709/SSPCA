# This function converts pconn data matrix into a vector
# @param pconn_data accepts the pconn data matrix
# @param pconn_vector accepts the pconn data in vector form
# @param direction coded to convert matrix into vector
# @keywords
# @export
# @examples
# Matrix2Vector(pconn_data=NULL,pconn_vector=NULL,direction="to_vector")

# Function pconn data into either matrix or vector form
Matrix2Vector <- function(pconn_data=NULL,pconn_vector=NULL,direction="to_vector") {
  if(direction=="to_matrix"){
    data_out <- pconn_data
    data_out[upper.tri(data_out)] <- pconn_vector
    data_out[lower.tri(data_out)] <- t(data_out)[lower.tri(data_out)]
  } else
  {
    data_out <- pconn_data[upper.tri(pconn_data)]
  }
  return(data_out)
}