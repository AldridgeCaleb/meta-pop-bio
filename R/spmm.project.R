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
#' block diagonal Movement matrix MM (see `spmm.project.matrix` for more details).
#' @param n_timesteps The number of time steps into the future that should be
#' projected.
#' @param n_stages The number of stages (rows) in the metapopulation state matrix
#' N.
#' @param n_patches The number of patches (columns) in the metapopulation state
#' matrix N.
#' @param ddf Density-dependent function parameters (see `?spmm.ddf.params`)
#' @param harv Additive (harvest) mortality. Can be element (applies to all) or a 
#' vector of added mortality equal in length to demographic matrix of row × 
#' column. 
#' @param deter A list of three vectors. The first two, `from` and `to`, identify
#' where deterrence of movement is made. The third, `d`, contains the proportions
#' by which movement is deterred. Currently deterrence is assumed equal for all
#' stages.
#' @param P vec-permutation matrix -- required if ddf, harv, or deter given
#' @param BB block diagonal demographic matrix -- required if ddf, harv, or deter given
#' @param MM block diagonal movement (dispersal) matrix -- required if ddf, harv, 
#' or deter given 
#'
#' @note
#' Ensure that the structural lh_orders of population vector `n` and projection
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
#' # Peregrine falcon example from Hunter and Caswell (2005), data from Wootton
#' # and Bell (1992). Continues example from `spmm.project.matrix`.
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
#' # Movement parameter values
#' dx1 <- 0.27  # only juveniles disperse
#' dx2 <- 1 - dx1
#' # Movement matrices for stages
#' Mx1 <- matrix(c(dx2, dx1, dx1, dx2), nrow = n_patches, byrow = TRUE)
#' Mx2 <- diag(x = 1, nrow = n_patches, ncol = n_patches)  # no Movement by adults
#' # Movement block matrix construction
#' MM <- blk.diag(list(Mx1, Mx2))
#'
#' # Arrangement by patches
#' group_by <- "patches"
#' # Assumed movement before demography
#' lh_order <- "move"
#'
#' # Projection matrix construction
#' A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
#' group_by = group_by, lh_order = lh_order)  # BB %*% t(P) %*% MM %*% P
#'
#' # Initial stages within patches (patch group_by)
#' n <- c(
#'   50, 22,  # Northern patch juveniles (row 1, column 1) and adults (row 2, column 1)
#'   40, 17   # Southern patch juveniles (row 1, column 2) and adults (row 2, column 2)
#' )
#' comment(n) <- "patches"  # vec comment attr for group_by
#'
#' # Number of time steps to project into the future
#' n_timesteps <- 50
#'
#' # Project spatial matrix population model
#' projs <- spmm.project(n, A, n_timesteps, n_stages, n_patches)
#'
#' @export
spmm.project <-
  function(n, A, n_timesteps,
           n_stages, n_patches, 
           ddf = NA, harv = NA, deter = NA,
           P, BB, MM) {
    
    try(if (is.null(comment(n)))
      stop("Please specify structure of n as either patches or stages (e.g., comment(n) <- 'patches'.')")
    )
    
    try(if (comment(n) != group_by)
      stop("Structure of n and A are not the same; both should include either 'patches' or 'stages'.")
    )
    
    lh_order <- comment(A)
    group_by <- strsplit(lh_order, " +")[[1]][1]
    A_lh_order <- strsplit(lh_order, " +")[[1]][2]
    
    if (lh_order == "patches demo") {
      mat <- matrix(nrow = n_stages * n_patches, ncol = n_timesteps)
      tryCatch({
        mat[, 1] <- as.vector(n)
      }, error = function(e) {
        stop("Length of n and n_stages × n_patches are not equal.")
      })
      
      if (any(!is.na(harv))) {
        matlist <- unblk.diag(BB, n_stages)
        for (i in seq_along(matlist)) {
          B <- matlist[i]
          M <-
            -log(B[[1]][-1, ])  # Transform survival probabilities to mortality rates
          if (is.vector(harv) & length(harv) == 1) {
            M <- M + harv  # Add constant harv to all mortality rates
          } else if (is.vector(harv) & length(harv) == length(matlist)){
            M <- M + harv[i]
          } else if (is.list(harv) & dim(harv[[i]]) == dim(matlist[[i]])) {
            M <- M + harv[[i]]  # Add list of harv mortality rates matrix-wise
          } else {
            stop(
              "Harvest mortality must be: a vector length == 1 if universally applied, a vector length == n_patches or n_stages (corresponding to group_by) if differential across grouping, or a list of matrices == length(BB) (corresponding to group_by) with dim(matrix) == dim(BB[[]]) if differential across grouping and in matrix elements."
            )
          }
          B[[1]][-1, ] <- exp(-M)  # Transform back to survival probabilities
          matlist[i] <- B
        }
        BB <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 group_by = group_by, lh_order = A_lh_order)
      }
      
      if (!is.na(deter)) {
        matlist <- unblk.diag(MM, n_patches)
        for (i in seq_along(matlist)) {
          M <- matlist[i]
          if (!is.identity.matrix(M)) {
            M[deter$from, deter$to] <- M[deter$from, deter$to] * deter$d
          }
          matlist[i] <- M
        }
        MM <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 group_by = group_by, lh_order = A_lh_order)
      }
      
      for (t in 2:n_timesteps) {
        if (!is.na(ddf)){
          matlist <- unblk.diag(BB, n_stages)
          for (i in seq_along(matlist)) {
            B <- matlist[i]
            if (ddf$f_type == "Ricker") {
              B[[1]][1, ] <- dd.rec.Ricker(mat[, t - 1], ddf$a[i], ddf$b[i], theta)
            } 
            if (ddf$f_type == "Beverton-Holt") {
              B[[1]][1, ] <- dd.rec.BevertonHolt(mat[, t - 1], ddf$a[i], ddf$b[i], theta)
            }
            if (ddf$s_type == "logistic") {
              B[[1]][1, ] <- dd.surv.logistic(mat[, t - 1], ddf$r[i], ddf$K[i])
            }
            if (ddf$s_type == "ddExponential") {
              B[[1]][1, ] <- dd.surv.exponential(mat[, t - 1], ddf$r[i], ddf$K[i])
            }
          }
          BB <- blk.diag(matlist)
          A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                   group_by = group_by, lh_order = A_lh_order)
        }

        # Projection to next t
        if (all(mat[, t - 1]%%1==0)) {
          mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
        } else {
          mat[, t] <- as.vector(A %*% mat[, t - 1]) 
        }
      }
      if (!is.null(rownames(n))) {
        rownames(mat) <- rownames(n)
      }
      colnames(mat) <- paste(1:n_timesteps)
      
    } else if (lh_order == "patches move") {
      mat <- matrix(nrow = n_stages * n_patches, ncol = n_timesteps)
      tryCatch({
        mat[, 1] <- as.vector(n)
      }, error = function(e) {
        stop("Length of n and n_stages × n_patches are not equal.")
      })
      
      if (any(!is.na(harv))) {
        matlist <- unblk.diag(BB, n_stages)
        for (i in seq_along(matlist)) {
          B <- matlist[i]
          M <-
            -log(B[[1]][-1, ])  # Transform survival probabilities to mortality rates
          if (is.vector(harv) & length(harv) == 1) {
            M <- M + harv  # Add constant harv to all mortality rates
          } else if (is.vector(harv) & length(harv) == length(matlist)){
            M <- M + harv[i]
          } else if (is.list(harv) & dim(harv[[i]]) == dim(matlist[[i]])) {
            M <- M + harv[[i]]  # Add list of harv mortality rates matrix-wise
          } else {
            stop(
              "Harvest mortality must be: a vector length == 1 if universally applied, a vector length == n_patches or n_stages (corresponding to group_by) if differential across grouping, or a list of matrices == length(BB) (corresponding to group_by) with dim(matrix) == dim(BB[[]]) if differential across grouping and in matrix elements."
            )
          }
          B[[1]][-1, ] <- exp(-M)  # Transform back to survival probabilities
          matlist[i] <- B
        }
        BB <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                 group_by = group_by, lh_order = A_lh_order)
      }
      
      if (!is.na(deter)) {
        matlist <- unblk.diag(MM, n_patches)
        for (i in seq_along(matlist)) {
          M <- matlist[i]
          if (!is.identity.matrix(M)) {
            M[deter$from, deter$to] <- M[deter$from, deter$to] * deter$d
          }
          matlist[i] <- M
        }
        MM <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 group_by = group_by, lh_order = A_lh_order)
      }
      
      for (t in 2:n_timesteps) {
        
        if (!is.na(ddf)){
          matlist <- unblk.diag(BB, n_stages)
          for (i in seq_along(matlist)) {
            B <- matlist[i]
            if (ddf$f_type == "Ricker") {
              a <- B[[1]][1, ]
              B[[1]][1, ] <- dd.rec.Ricker(mat[, t - 1], a, b)
            } 
            if (ddf$f_type == "Beverton-Holt") {
              B[[1]][1, ] <- B[[1]][1, ] * dd.rec.BevertonHolt(mat[, t - 1], ddf$r[i], ddf$K[i])
            }
            if (ddf$s_type == "logistic") {
              B[[1]][-1, ] <- B[[1]][-1, ] * dd.surv.logistic(mat[, t - 1], ddf$r[i], ddf$K[i])
            }
            matlist[i] <- B
          }
          BB <- blk.diag(matlist)
          A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                   group_by = group_by, lh_order = A_lh_order)
        }
        
        # Projection to next t
        if (all(mat[, t - 1]%%1==0)) {
          mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
        } else {
          mat[, t] <- as.vector(A %*% mat[, t - 1]) 
        }
      }
      if (!is.null(rownames(n))) {
        rownames(mat) <- rownames(n)
      }
      colnames(mat) <- paste(1:n_timesteps)
      
    } else if (lh_order == "stages demo") {
      mat <- matrix(nrow = n_patches * n_stages, ncol = n_timesteps)
      tryCatch({
        mat[, 1] <- as.vector(n)
      }, error = function(e) {
        stop("Length of n and n_stages × n_patches are not equal.")
      })
      
      if (any(!is.na(harv))) {
        matlist <- unblk.diag(BB, n_stages)
        for (i in seq_along(matlist)) {
          B <- matlist[i]
          M <-
            -log(B[[1]][-1, ])  # Transform survival probabilities to mortality rates
          if (is.vector(harv) & length(harv) == 1) {
            M <- M + harv  # Add constant harv to all mortality rates
          } else if (is.vector(harv) & length(harv) == length(matlist)){
            M <- M + harv[i]
          } else if (is.list(harv) & dim(harv[[i]]) == dim(matlist[[i]])) {
            M <- M + harv[[i]]  # Add list of harv mortality rates matrix-wise
          } else {
            stop(
              "Harvest mortality must be: a vector length == 1 if universally applied, a vector length == n_patches or n_stages (corresponding to group_by) if differential across grouping, or a list of matrices == length(BB) (corresponding to group_by) with dim(matrix) == dim(BB[[]]) if differential across grouping and in matrix elements."
            )
          }
          B[[1]][-1, ] <- exp(-M)  # Transform back to survival probabilities
          matlist[i] <- B
        }
        BB <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 group_by = group_by, lh_order = A_lh_order)
      }
      
      if (!is.na(deter)) {
        matlist <- unblk.diag(MM, n_patches)
        for (i in seq_along(matlist)) {
          M <- matlist[i]
          if (!is.identity.matrix(M)) {
            M[deter$from, deter$to] <- M[deter$from, deter$to] * deter$d
          }
          matlist[i] <- M
        }
        MM <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                 group_by = group_by, lh_order = A_lh_order)
      }
      
      for (t in 2:n_timesteps) {
        if (!is.na(ddf)){
          matlist <- unblk.diag(BB, n_stages)
          for (i in seq_along(matlist)) {
            B <- matlist[i]
            if (ddf$f_type == "Ricker") {
              B[[1]][1, ] <- B[[1]][1, ] * dd.rec.Ricker(mat[, t - 1], ddf$r[i], ddf$K[i])
            } 
            if (ddf$f_type == "Beverton-Holt") {
              B[[1]][1, ] <- B[[1]][1, ] * dd.rec.BevertonHolt(mat[, t - 1], ddf$r[i], ddf$K[i])
            }
            if (ddf$s_type == "logistic") {
              B[[1]][-1, ] <- B[[1]][-1, ] * dd.surv.logistic(mat[, t - 1], ddf$r[i], ddf$K[i])
            }
            matlist[i] <- B
          }
          BB <- blk.diag(matlist)
          A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                   group_by = group_by, lh_order = A_lh_order)
        }
        
        # Projection to next t
        if (all(mat[, t - 1]%%1==0)) {
          mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
        } else {
          mat[, t] <- as.vector(A %*% mat[, t - 1]) 
        }
      }
      if (!is.null(rownames(n))) {
        rownames(mat) <- rownames(n)
      }
      colnames(mat) <- paste(1:n_timesteps)
      
    } else if (lh_order == "stages move") {
      mat <- matrix(nrow = n_patches * n_stages, ncol = n_timesteps)
      tryCatch({
        mat[, 1] <- as.vector(n)
      }, error = function(e) {
        stop("Length of n and n_stages × n_patches are not equal.")
      })
      
      if (any(!is.na(harv))) {
        matlist <- unblk.diag(BB, n_stages)
        for (i in seq_along(matlist)) {
          B <- matlist[i]
          M <-
            -log(B[[1]][-1, ])  # Transform survival probabilities to mortality rates
          if (is.vector(harv) & length(harv) == 1) {
            M <- M + harv  # Add constant harv to all mortality rates
          } else if (is.vector(harv) & length(harv) == length(matlist)){
            M <- M + harv[i]
          } else if (is.list(harv) & dim(harv[[i]]) == dim(matlist[[i]])) {
            M <- M + harv[[i]]  # Add list of harv mortality rates matrix-wise
          } else {
            stop(
              "Harvest mortality must be: a vector length == 1 if universally applied, a vector length == n_patches or n_stages (corresponding to group_by) if differential across grouping, or a list of matrices == length(BB) (corresponding to group_by) with dim(matrix) == dim(BB[[]]) if differential across grouping and in matrix elements."
            )
          }
          B[[1]][-1, ] <- exp(-M)  # Transform back to survival probabilities
          matlist[i] <- B
        }
        BB <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                 group_by = group_by, lh_order = A_lh_order)
      }
      
      if (!is.na(deter)) {
        matlist <- unblk.diag(MM, n_patches)
        for (i in seq_along(matlist)) {
          M <- matlist[i]
          if (!is.identity.matrix(M)) {
            M[deter$from, deter$to] <- M[deter$from, deter$to] * deter$d
          }
          matlist[i] <- M
        }
        MM <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                 group_by = group_by, lh_order = A_lh_order)
      }
      
      for (t in 2:n_timesteps) {
        if (!is.na(ddf)){
          matlist <- unblk.diag(BB, n_stages)
          for (i in seq_along(matlist)) {
            B <- matlist[i]
            if (ddf$f_type == "Ricker") {
              B[[1]][1, ] <- B[[1]][1, ] * dd.rec.Ricker(mat[, t - 1], ddf$r[i], ddf$K[i])
            } 
            if (ddf$f_type == "Beverton-Holt") {
              B[[1]][1, ] <- B[[1]][1, ] * dd.rec.BevertonHolt(mat[, t - 1], ddf$r[i], ddf$K[i])
            }
            if (ddf$s_type == "logistic") {
              B[[1]][-1, ] <- B[[1]][-1, ] * dd.surv.logistic(mat[, t - 1], ddf$r[i], ddf$K[i])
            }
            matlist[i] <- B
          }
          BB <- blk.diag(matlist)
          A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                   group_by = group_by, lh_order = A_lh_order)
        }
        
        # Projection to next t
        if (all(mat[, t - 1]%%1==0)) {
          mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
        } else {
          mat[, t] <- as.vector(A %*% mat[, t - 1]) 
        }
      }
      if (!is.null(rownames(n))) {
        rownames(mat) <- rownames(n)
      }
      colnames(mat) <- paste(1:n_timesteps)
      
    }
    
    if (A_lh_order == "move") {
      A_TYPE <- "movement then demography"
    } else if (A_lh_order == "demo") {
      A_TYPE <- "demography then movement"
    }
    print(
      paste(
        "Deterministic spatial matrix model projections for",
        group_by,
        "structured population vector and",
        A_TYPE,
        "A projection matrix."
      )
    )
    if (group_by == "patches") {
      comment(mat) <- paste(c(
        group_by,
        as.character(n_stages),
        as.character(n_patches)
      ))
    } else if (group_by == "stages") {
      comment(mat) <- paste(c(
        group_by,
        as.character(n_patches),
        as.character(n_stages)
      ))
    }
    return(mat)
  }
