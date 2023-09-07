#' @title vec-permutation matrix
#' 
#' @description
#' Stacks the columns of a matrix one on top of the other to covert a population 
#' vector organized by patches into a population vector organized by stages and 
#' vice versa.
#' 
#' For use in spatial matrix population models described by Hunter and Caswell
#' (2005) and Lebreton (1996).
#' 
#' @param n_stages The number of stages (rows) in the metapopulation state matrix 
#' N.
#' @param n_patches The number of patches (columns) in the metapopulation state 
#' matrix N.
#' @param group_by The structural form of the population vector `n` as either
#' "patches" or "stages".
#' 
#' @note
#' The population vector n can be written in two ways: by patches or by stages. 
#' If modeling demography then dispersal, it is convenient to construct n by 
#' patches. If modeling dispersal then demography, it is convenient to construct 
#' n by stages. Please refer to Hunter and Caswell (2005) for details.
#' 
#' For convenience, the argument `group_by` is used to properly construct the
#' vec-permutation matrix correctly.
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
#' and Bell (1992).
#' 
#' Define the number of patches and stages
#' n_patches <- 2  # northern = 1x; southern = 2x
#' n_stages <- 2  # juvenile = x1; adult = x2
#' group_by <- "patches"
#' 
#' Construct vec-permutation matrix
#' P <- vec.perm(n_stages, n_patches, group_by)
#' 
#' @export
vec.perm <-
  function(stages,
           patches,
           group_by = c("patches", "stages")) {
    if (group_by == "patches") {
      m <- stages
      n <- patches
      P <- matrix(0, nrow = m * n, ncol = m * n)
      for (i in 1:m) {
        for (j in 1:n) {
          row_idx <- (i - 1) * n + j
          col_idx <- (j - 1) * m + i
          P[row_idx, col_idx] <- 1
        }
      }
    } else if (group_by == "stages") {
      m <- patches
      n <- stages
      P <- matrix(0, nrow = m * n, ncol = m * n)
      for (i in 1:m) {
        for (j in 1:n) {
          row_idx <- (i - 1) * n + j
          col_idx <- (j - 1) * m + i
          P[row_idx, col_idx] <- 1
        }
      }
    }
    comment(P) <- group_by
    return(P)
  }
