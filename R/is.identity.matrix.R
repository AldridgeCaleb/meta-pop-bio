#' @title TRUE/FALSE test for identity matrix
#' 
#' @param mat A matrix.
#' 
#' @return TRUE/FALSE
#' 
#' @examples
#' M <- matrix(nrow = 3)
#' is.identity.matrix(M)
#' 
#' @export
is.identity.matrix <- function(mat) {
  if (nrow(mat) != ncol(mat)) {
    return(FALSE)
  }
  identity_mat <- diag(nrow = nrow(mat))
  return(identical(mat, identity_mat))
}
