#' @title sensitivity of projection matrix A
#' 
#' @description
#' Calculates the sensitivity of lambda to changes in the projection matrix A 
#' (see `spmm.proj.matrix`) as described in Hunter and Caswell (2005) and Lebreton (1996). 
#' 
#' @param A The spatial population projection matrix constructed from the
#' vec-permutation matrix P, block diagonal demographic matrix BB, and 
#' block diagonal movement matrix MM (see `spmm.proj.matrix` for more details).
#' 
#' @returns A matrix containing sensitivity values for the projection matrix A.
#' According to Morris and Doak (2003) sensitivity values od lambda for a 
#' particular matrix element is "directly proportional to the fraction of 
#' individuals in the population on which the element will act times the future 
#' value of each individual that the element 'creates'" (p. 226). 
#'  
#' @note
#' Ensure that the structural types of population vector `n` and projection 
#' matrix `A` are the same. Otherwise, projections may produce incorrect values!
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
#' Peregrine falcon example from Hunter and Caswell (2005), data from Wootton
#' and Bell (1992). Continues example from `spmm.proj.matrix`.
#' 
#' Define the number of patches and stages
#' n_patches <- 2  # northern = 1x; southern = 2x
#' n_stages <- 2  # juvenile = x1; adult = x2
#' group_by <- "patches"
#' 
#' Construct vec-permutation matrix
#' P <- vec.perm(n_stages, n_patches, group_by)
#' 
#' Demographic parameter values
#' Northern
#' f11 <- 0.00  # only adults reproduce
#' f12 <- 0.26
#' s11 <- 0.72
#' s12 <- 0.77
#' Southern
#' f21 <- 0.00
#' f22 <- 0.19  
#' s21 <- 0.72
#' s22 <- 0.77
#' 
#' Demography matrices for patches
#' B1x <-
#'   matrix(c(f11, f12, s11, s12),
#'          nrow = 2,
#'          byrow = TRUE)
#' B2x <-
#'   matrix(c(f21, f22, s21, s22),
#'          nrow = 2,
#'          byrow = TRUE)
#' Demography block matrix construction
#' BB <- blk.diag(B1x, B2x)
#' 
#' Movement parameter values
#' dx1 <- 0.27  # only juveniles disperse
#' dx2 <- 1 - dx1
#' Movement matrices for stages
#' Mx1 <- matrix(c(dx2, dx1, dx1, dx2), nrow = n_patches, byrow = TRUE)
#' Mx2 <- diag(x = 1, nrow = n_patches, ncol = n_patches)  # no movement by adults
#' Movement block matrix construction
#' MM <- blk.diag(Mx1, Mx2)
#' 
#' Arrangement by patches
#' group_by <- "patches"
#' Assumed movement before demography
#' type <- "move"
#' 
#' Projection matrix construction
#' A <- spmm.proj.matrix(P, BB, MM, group_by, type)  # BB %*% t(P) %*% MM %*% P 
#' 
#' Calculate sensitivity of lambda to elements of projection matrix A
#' A_sens <- spmm.proj.matrix.sens(A)
#' 
#' @export
spmm.proj.matrix.sens <- function(A) {
  eig <- eigen(A)
  lambda <- max(Re(eig$values))
  v <- eig$vectors[, which.max(Re(eig$values))]
  eig_t <- eigen(t(A))
  w <- eig_t$vectors[, which.max(Re(eig_t$values))]
  SA <- t((v %*% t(w)) / sum(w * v))
  return(SA)
}
