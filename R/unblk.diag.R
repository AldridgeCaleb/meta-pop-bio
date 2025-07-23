#' @title unblock a diagonal matrix
#' 
#' @description
#' Deconstructs a block matrix into component matrices. 
#' 
#' @param blk_matrix List of matrices, separated by commas, for constructing the 
#' block diagonal matrix.
#' @param dimensions Number of stages (n_stages) or patches (n_patches), depending on if a demographic 
#' or movement (dispersal) block matrix is being decomposed.
#' 
#' @examples
#' # Peregrine falcon example from Hunter and Caswell (2005), data from Wootton
#' # and Bell (1992). Continues example from
#' # `vec.perm`.
#' 
#' # Define the number of patches and stages
#' n_patches <- 2  # northern = 1x; southern = 2x
#' n_stages <- 2  # juvenile = x1; adult = x2
#' group_by <- "patches"
#' 
#' # Construct vec-permutation matrix
#' P <- vec.perm(n_stages, n_patches, group_by)
#' 
#' # Demographic parameter values
#' # Northern
#' f11 <- 0.00  # only adults reproduce
#' f12 <- 0.26
#' s11 <- 0.72
#' s12 <- 0.77
#' # Southern
#' f21 <- 0.00
#' f22 <- 0.19  
#' s21 <- 0.72
#' s22 <- 0.77
#' 
#' # Demography matrices for patches
#' B1x <-
#'   matrix(c(f11, f12, s11, s12),
#'          nrow = 2,
#'          byrow = TRUE)
#' B2x <-
#'   matrix(c(f21, f22, s21, s22),
#'          nrow = 2,
#'          byrow = TRUE)
#' # Demography block matrix construction
#' BB <- blk.diag(list(B1x, B2x))
#' 
#' # Dispersal parameter values
#' dx1 <- 0.27  # only juveniles disperse
#' dx2 <- 1 - dx1
#' # Dispersal matrices for stages
#' Mx1 <- matrix(c(dx2, dx1, dx1, dx2), nrow = n_patches, byrow = TRUE)
#' Mx2 <- diag(x = 1, nrow = n_patches, ncol = n_patches)  # no dispersal by adults
#' # Dispersal block matrix construction
#' MM <- blk.diag(list(Mx1, Mx2))
#' 
#' matlist <- unblk.diag(MM, n_patches)
#' 
#' @export
unblk.diag  <- function(blk_matrix, dimensions) {
  n_mats <- ncol(blk_matrix) / dimensions 
  
  matlist <- list()

  # Loop over each block and extract the submatrix
  for (i in seq_len(n_mats)) {
    start_idx <- (i - 1) * dimensions + 1
    end_idx <- start_idx + dimensions - 1
    sub_mat <- blk_matrix[start_idx:end_idx, 
                          start_idx:end_idx]
    matlist[[i]] <- sub_mat
  }
  
  return(matlist)
}
