#' @title projection matrix
#' 
#' @description
#' Constructs projection matrix A given a vec-permutation matrix P, a 
#' demographic block diagonal matrix BB, and a movement (dispersal) block diagonal
#' matrix MM as described in Hunter and Caswell (2005) and Lebreton (1996). 
#' 
#' @param P The vec-permutation matrix (see `vec.perm` function). 
#' @param BB The block diagonal demographics matrix (see `blk.diag` function).
#' @param MM The block diagonal movement matrix (see `blk.diag` function). 
#' @param group_by The structural form of the population vector `n` as either
#' "patches" or "stages".
#' @param lh_order The order of demographic and movement (dispersal) processes.
#' 
#' @note
#' The demographic block matrix BB consists of s × s number of 
#' demographic matrices (one for each patch). The demographic block matrix MM 
#' consists of p × p number of stage movement matrices (one for each 
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
#' # and Bell (1992). Continues example from `blk.diag`.
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
#' # Movement matrices for patches
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
#' # Movement parameter values
#' dx1 <- 0.27  # only juveniles disperse
#' dx2 <- 1 - dx1
#' # Movement matrices for stages
#' Mx1 <- matrix(c(dx2, dx1, dx1, dx2), nrow = n_patches, byrow = TRUE)
#' Mx2 <- diag(x = 1, nrow = n_patches, ncol = n_patches)  # no movement by adults
#' # Movement block matrix construction
#' MM <- blk.diag(list(Mx1, Mx2))
#' 
#' # Arrangement by patches
#' group_by <- "patches"
#' # Assumed movement before demography
#' lh_order <- "move"
#' 
#' # Projection matrix construction
#' A <- spmm.project.matrix(P, BB, MM, group_by, lh_order)  # BB %*% t(P) %*% MM %*% P 
#' 
#' @export
spmm.project.matrix <- function(P, BB, MM, group_by = c("patches", "stages"), 
                                lh_order = c("demo", "move")) {
  A <- NULL  # Initialize A to NULL
  if (group_by == "patches" && lh_order == "demo") {
    A <- t(P) %*% MM %*% P %*% BB 
  } else if (group_by == "patches" && lh_order == "move") {
    A <- BB %*% t(P) %*% MM %*% P
  } else if (group_by == "stages" && lh_order == "demo") {
    A <- MM %*% P %*% BB %*% t(P)
  } else if (group_by == "stages" && lh_order == "move") {
    A <- P %*% BB %*% t(P) %*% MM
  } else {
    stop("Invalid combination of group_by and lh_order.")  # Default case
  }
  comment(A) <- paste(group_by, lh_order)
  return(A)
}
