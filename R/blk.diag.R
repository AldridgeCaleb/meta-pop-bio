#' @title block diagonal matrix construction
#' 
#' @description
#' Constructs a block matrix of of matrices. Matrices types for spatial matrix
#' population models described in Hunter and Caswell (2005) and Lebreton (1996), 
#' include demographic B and dispersal M for which block diagonal matrices are 
#' constructed, BB and MM for demographic and dispersal respectively.
#' 
#' @param matlist List of matrices, separated by commas, for constructing the 
#' block diagonal matrix.
#' 
#' @note
#' The demographic block matrix BB consists of s × s number of 
#' demographic matrices (one for each patch). The demographic block matrix MM 
#' consists of p × p number of stage dispersal matrices (one for each 
#' stage). Please refer to Hunter and Caswell (2005) for details.
#' 
#' @references
#' Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988). The New S Language. 
#' Wadsworth & Brooks/Cole.
#' 
#' Wootton, J.T., and Bell, D.A. (1992). A metapopulation model of the peregrine 
#' falcon in California: viability and management strategies. Ecological 
#' Applications 2:307--321.
#' 
#' Lebreton, J. D. (1996). Demographic models for subdivided populations: the 
#' renewal equation approach. Theoretical Population Biology 49:291--313.
#' 
#' Caswell, H. (2001). Matrix Population Models: Construction, analysis, and 
#' interpretation (2nd ed.). Sinauer Associates.
#' 
#' Morris, W. F., and Doak, D. F. (2003). Quantitative Conservation Biology: 
#' Theory and practice of population viability analysis. Sinauer Associates.
#' 
#' Hunter, C. M. and Caswell, H. (2005). The use of vec-permutation matrix in
#' spatial matrix population models. Ecological Modelling 188:15--21.
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
#' @export
blk.diag  <- function(matlist) {
  if (!is.list(matlist)) {
    Ms <- list(...)
  } else if (is.null(matlist)) {
    stop("Please provide list of matrices!")
  }
  Ms <- matlist
  n_rows <- sum(sapply(Ms, nrow))
  n_cols <- sum(sapply(Ms, ncol))
  MM <- matrix(0, n_rows, n_cols)
  m <- 1
  n <- 1
  for (i in Ms) {
    mm <- m + nrow(i) - 1
    nn <- n + ncol(i) - 1
    MM[m:mm, n:nn] <- i
    m <- mm + 1
    n <- nn + 1
  }
  return(MM)
}
