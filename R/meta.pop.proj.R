#' @title deterministic projection of spatial matrix population model
#' 
#' @description
#' Produces deterministic future population sizes of stages s in patches p.
#' Projections can be ordered with demography then movement (dispersal) or vice 
#' versa using population vectors arranged by patches or stages as described in 
#' Hunter and Caswell (2005) and Lebreton (1996). 
#' 
#' @param n Vector of stage and patch starting values. If structure for `n` is
#' by patches, the vector sould order values for each stage within each patch 
#' and "stack" patches. If structure for `n` is by stages, the vector sould 
#' order values for each patch within each stage and "stack" stages.
#' @param A The spatial population projection matrix constructed from the
#' vec-permutation matrix P, block diagonal demographic matrix BB, and 
#' block diagonal Movement matrix MM (see `proj.mat` for more details).
#' @param n_timesteps The number of time steps into the future that should be 
#' projected.
#' @param n_stages The number of stages (rows) in the metapopulation state matrix 
#' N.
#' @param n_patches The number of patches (columns) in the metapopulation state 
#' matrix N.
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
#' Applications 2:307-–321.
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
#' and Bell (1992). Continues example from `proj.mat`.
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
#' Mx2 <- diag(x = 1, nrow = n_patches, ncol = n_patches)  # no Movement by adults
#' Movement block matrix construction
#' MM <- blk.diag(Mx1, Mx2)
#' 
#' Arrangement by patches
#' group_by <- "patches"
#' Assumed movement before demography
#' type <- "move"
#' 
#' Projection matrix construction
#' A <- proj.mat(P, BB, MM, group_by, type)  # BB %*% t(P) %*% MM %*% P 
#' 
#' Initial stages within patches (patch group_by)  
#' n <- c(
#'   50, 22,  # Northern patch adults and juveniles
#'   40, 17   # Southern patch adults and juveniles
#' )
#' comment(n) <- "patches"  # vec comment attr for group_by 
#' 
#' Number of time steps to project into the future
#' n_timesteps <- 50
#' 
#' Project spatial matrix population model
#' projs <- meta.pop.proj(n, A, n_timesteps, n_stages, n_patches)
#' 
#' @export
meta.pop.proj <- function(n, A, n_timesteps, n_stages, n_patches) {
  type <- comment(A)
  group_by <- strsplit(type, " +")[[1]][1]
  A_type <- strsplit(type, " +")[[1]][2]
  
  try(if (is.null(comment(n)))
    stop("Please specify structure of n as either patches or stages (e.g., comment(n) <- 'patches'.')")
  )
  
  try(if (comment(n) != group_by)
    stop("Structure of n and A are not the same; both should include either 'patches' or 'stages'.")
  )
    
  if (type == "patches demo") {
    mat <- matrix(nrow = n_stages * n_patches, ncol = n_timesteps)
    # throw.err <-
    #   function(...)
    #     stop("Length of n and n_stages × n_patches are not equal.")
    # try(throw.err(
      mat[, 1] <- as.vector(n)
      # ))
    
    for (t in 2:n_timesteps){
      mat[, t] <- as.vector(A %*% mat[, t - 1])
    }
    if (!is.null(rownames(n))) {
      rownames(mat) <- rownames(n)
    }
    colnames(mat) <- paste(1:n_timesteps)
    
  } else if (type == "patches move") {
    mat <- matrix(nrow = n_stages * n_patches, ncol = n_timesteps)
    tryCatch({
      mat[, 1] <- as.vector(n)
    }, error = function(e) {
      stop("Length of n and n_stages × n_patches are not equal.")
    })
    
    for (t in 2:n_timesteps){
      mat[, t] <- as.vector(A %*% mat[, t - 1])
    }
    if (!is.null(rownames(n))) {
      rownames(mat) <- rownames(n)
    }
    colnames(mat) <- paste(1:n_timesteps)
    
  } else if (type == "stages demo") {
    mat <- matrix(nrow = n_patches * n_stages, ncol = n_timesteps)
    tryCatch({
      mat[, 1] <- as.vector(n)
    }, error = function(e) {
      stop("Length of n and n_stages × n_patches are not equal.")
    })
    
    for (t in 2:n_timesteps){
      mat[, t] <- as.vector(A %*% mat[, t - 1])
    }
    if (!is.null(rownames(n))) {
      rownames(mat) <- rownames(n)
    }
    colnames(mat) <- paste(1:n_timesteps)
    
  } else if (type == "stages move") {
    mat <- matrix(nrow = n_patches * n_stages, ncol = n_timesteps)
    # throw.err <-
    #   function(...)
    #     stop("Length of n and n_stages × n_patches are not equal.")
    # try(throw.err(
      mat[, 1] <- as.vector(n)
      # ))
    
    for (t in 2:n_timesteps){
      mat[, t] <- as.vector(A %*% mat[, t - 1])
    }
    if (!is.null(rownames(n))) {
      rownames(mat) <- rownames(n)
    }
    colnames(mat) <- paste(1:n_timesteps)
    
  }
  
  if (A_type == "move") {
    A_TYPE <- "movement then demography"
  } else if (A_type == "demo") {
    A_TYPE <- "demography then movement"
  }
  print(paste("Deterministic spatial matrix model projections for", group_by, 
              "structured population vector and", A_TYPE, "A projection matrix."))
  if (group_by == "patches") {
    comment(mat) <- paste(c(group_by,
                            as.character(n_stages), 
                            as.character(n_patches)))
  } else if (group_by == "stages") {
    comment(mat) <- paste(c(group_by,
                            as.character(n_patches), 
                            as.character(n_stages)))
  }
  return(mat)
} 
