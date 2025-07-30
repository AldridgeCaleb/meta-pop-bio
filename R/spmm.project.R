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
#' @param mod_mort Optional modification to survival rates via additive mortality. 
#' Can be used to increase mortality (decrease survival) across one or more 
#' demographic matrices. This argument accepts one of the following forms:
#'   - **Single numeric value**: A scalar applied uniformly to all non-recruitment 
#'   elements across all demographic matrices (i.e., same additive mortality everywhere).
#'   - **Numeric vector**: A vector of length equal to either `n_patches` or `n_stages`, 
#'   depending on `group_by`. Each value is applied to the corresponding demographic 
#'   matrix in the list.
#'   - **List of matrices**: A list of the same length as the number of demographic 
#'   matrices (`n_patches` or `n_stages`), with each matrix matching the dimensions 
#'   of its corresponding demographic matrix (excluding the recruitment row). This 
#'   allows specifying element-wise additive mortality for each matrix.
#' Use a single value for uniform mortality, a vector for matrix-wise variation, 
#' or a list of matrices for full element-wise control.
#' @param mod_rec Optional modification to recruitment rates via additive effects. 
#' This argument can be specified in one of three forms, similar to `mod_mort`:
#'   - **Single numeric value**: A scalar applied uniformly to the recruitment row 
#'   (top row) of all demographic matrices.
#'   - **Numeric vector**: A vector of length equal to either `n_patches` or `n_stages` 
#'   (depending on `group_by`), where each value modifies the recruitment row across 
#'   the entire corresponding demographic matrix.
#'   - **List of numeric vectors**: A list of length equal to the number of 
#'   demographic matrices (`n_patches` or `n_stages`), where each list element is 
#'   a vector of length equal to the number of columns (`ncol`) in the demographic 
#'   matrix. This allows element-wise modification of recruitment rates within each 
#'   matrix.
#' Use a single value for uniform recruitment modification, a vector for matrix-wise 
#' differences, or a list of vectors for full element-wise control across recruitment 
#' elements.
#' @param mod_move Optional modification to movement. This is specified as a data.frame
#' that must contain columns:
#'     - `stage`: the index of the stage to which the modification applies,
#'     - `triangle`: either `"lower"` or `"upper"`, specifying which half of the 
#'     movement matrix to modify,
#'     - `at`: the column index of the movement matrix where modification is applied,
#'     - `d_perc`: the proportional change (e.g., `0.5` for 50% reduction, `1.5` 
#'     for 50% increase) to apply to all off-diagonal elements in the specified 
#'     triangle and column.
#' After applying the modification, the corresponding column is re-normalized to 
#' ensure they sum to 1. Note: specific patch-to-patch modifications using `from`
#' and `to` columns is being developed for future functionality.
#' Currently modifications are currently applied equally across all stages.
#' @param P (optional) vec-permutation matrix -- required if ddf, mod_mort, or mod_move given
#' @param BB (optional) block diagonal demographic matrix -- required if ddf, mod_mort, or mod_move given
#' @param MM (optional) block diagonal movement (dispersal) matrix -- required if ddf, mod_mort, 
#' or mod_move given 
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
           ddf = NULL, mod_mort = NA, 
           mod_rec = NA, mod_move = NA,
           P, BB, MM) {

# Initial error handling --------------------------------------------------
    lh_order <- comment(A)
    group_by <- strsplit(lh_order, " +")[[1]][1]
    A_lh_order <- strsplit(lh_order, " +")[[1]][2]
    
    try(if (is.null(comment(n)))
      stop("Please specify structure of n as either patches or stages (e.g., comment(n) <- 'patches'.')")
    )
    
    try(if (comment(n) != group_by)
      stop("Structure of n and A are not the same; both should include either 'patches' or 'stages'.")
    )
    
# LH ORDER: patches - demo ------------------------------------------------
    if (lh_order == "patches demo") {
      mat <- matrix(nrow = n_stages * n_patches, ncol = n_timesteps)
      tryCatch({
        mat[, 1] <- as.vector(n)
      }, error = function(e) {
        stop("Length of n and n_stages × n_patches are not equal.")
      })
## Mortality modification
      if (any(!is.na(mod_mort))) {
        matlist <- unblk.diag(BB, n_stages)
        if (is.vector(mod_mort) &&
            length(mod_mort) > 1 && length(mod_mort) != length(matlist)) {
          stop("mod_mort vector must be length 1 or equal to the number of demographic matrices.")
        } 
        
        if (is.list(mod_mort)) {
          if (length(mod_mort) != length(matlist)) {
            stop("mod_mort list length must match number of demographic matrices.")
          } 
          
          for (j in seq_along(mod_mort)) {
            if (!all(dim(mod_mort[[j]]) == dim(matlist[[j]]) - c(1, 0))) {
              stop(paste("mod_mort matrix at index", j, "does not match dimensions of the demographic matrix."))
            }
          }
        }
        for (i in seq_along(matlist)) {
          B <- matlist[i]
          M <- -log(B[[1]][-1, ])  # Transform survival probabilities to mortality rates
          if (is.vector(mod_mort) && length(mod_mort) == 1) {
            M <- M + mod_mort
          } else if (!is.list(mod_mort) &&
                     length(mod_mort) == length(matlist)) {
            M <- M + mod_mort[i]
          } else if (is.list(mod_mort)) {
            M <- M + mod_mort[[i]]
          }
          B[[1]][-1, ] <- exp(-M)
          matlist[i] <- B
        }
        BB <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 group_by = group_by, 
                                 lh_order = A_lh_order)
      }
## Movement modification
      if (all(!is.na(mod_move)) && is.data.frame(mod_move)) {
        matlist <- unblk.diag(MM, n_patches)
        for (i in seq_along(matlist)) {
          M <- matlist[i]
          mm <- mod_move[mod_move$stage == i, ]
          for (j in unique(mm$triangle)) {
            if (j == "lower") {
              if ("at" %in% colnames(mm)) {
                for (k in unique(mm$at)) {
                  M[lower.tri(M) & (col(M) == k)] <- 
                    M[lower.tri(M) & (col(M) == k)] * mm[mm$at == k, "d_perc"]
                  M[k, k] <- M[k, k] + (1 - sum(M[, k]))
                }
              } else if (all(c("from", "to") %in% colnames(mm))) {
                stop("Funcitnality for from - to specification coming soon.")
              }
            } else if (j == "upper") {
              if ("at" %in% colnames(mm)) {
                for (k in unique(mm$at)) {
                  M[upper.tri(M) & (col(M) == k)] <- 
                    M[upper.tri(M) & (col(M) == k)] * mm[mm$at == k, "d_perc"]
                  M[k, k] <- M[k, k] + (1 - sum(M[, k]))
                }
              } else if (all(c("from", "to") %in% colnames(mm))) {
                stop("Funcitnality for from - to specification coming soon.")
              }
            }
          }
          matlist[i] <- M
        }
        MM <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 group_by = group_by,
                                 lh_order = A_lh_order)
      }
# ## Density-dependence
#       for (t in 2:n_timesteps) {
#         if (any(!is.null(ddf))){
#           matlist <- unblk.diag(BB, n_stages)
#           for (i in seq_along(matlist)) {
#             B <- matlist[i]
#             if (ddf$type == "Ricker") {
#               B[[1]][1, ] <- dd.rec.Ricker(mat[, t - 1], ddf$a[i], ddf$b[i], theta)
#             } 
#             if (ddf$type == "Beverton-Holt") {
#               B[[1]][1, ] <- dd.rec.BevertonHolt(mat[, t - 1], ddf$a[i], ddf$b[i], theta)
#             }
#             if (ddf$type == "logistic") {
#               B[[1]][1, ] <- dd.growth.logistic(N = mat[c(i * n_stages - 1):c(i * n_stages), t - 1], 
#                                                 r = ddf$r[i], B = B[[1]][1, ], K = ddf$K[i])
#             }
#             if (ddf$type == "ddExponential") {
#               B[[1]][1, ] <- dd.growth.exponential(mat[, t - 1], ddf$r[i], ddf$K[i])
#             }
#             if (ddf$type == "general") {
#               B[[1]][1, ] <- dd.growth.general(mat[, t - 1], ddf$r[i], ddf$K[i], ddf$theta)
#             }
#             matlist[i] <- B
#           }
#           BB <- blk.diag(matlist)
#           A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
#                                    group_by = group_by, lh_order = A_lh_order)
#         }
# ## Projection
#         if (all(mat[, t - 1]%%1==0)) {
#           mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
#         } else {
#           mat[, t] <- as.vector(A %*% mat[, t - 1]) 
#         }
#       }
#       if (!is.null(rownames(n))) {
#         rownames(mat) <- rownames(n)
#       }
#       colnames(mat) <- paste(1:n_timesteps)
      
      for (t in 2:n_timesteps) {
        if (!is.null(ddf)) {
          matlist <- unblk.diag(BB, n_stages)
          for (i in seq_along(matlist)) {
            B <- matlist[i]
            idx <- ((i - 1) * n_stages + 1):(i * n_stages)
            Ni <- mat[idx, t - 1]
            B[[1]][1, ] <- dd.growth.logistic(N = Ni, B = B[[1]][1, ], K = ddf$K[i],
                                              beta = ddf$beta, theta = ddf$theta)
            matlist[i] <- B
          }
          BB <- blk.diag(matlist)
          A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                   group_by = group_by, lh_order = A_lh_order)
        }
        if (anyNA(mat[, t - 1])) {
          warning(paste("NA detected at timestep", t - 1))
          print(which(is.na(mat[, t - 1])))
        }
        
        if (all(mat[, t - 1] %% 1 == 0)) {
          mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
        } else {
          mat[, t] <- as.vector(A %*% mat[, t - 1])
        }
      }
      if (!is.null(rownames(n))) {
        rownames(mat) <- rownames(n)
      }
      colnames(mat) <- paste(1:n_timesteps)
      
# LH ORDER: patches - move ------------------------------------------------
    } else if (lh_order == "patches move") {
      mat <- matrix(nrow = n_stages * n_patches, ncol = n_timesteps)
      tryCatch({
        mat[, 1] <- as.vector(n)
      }, error = function(e) {
        stop("Length of n and n_stages × n_patches are not equal.")
      })
## Mortality modification
      if (any(!is.na(mod_mort))) {
        matlist <- unblk.diag(BB, n_stages)
        if (is.vector(mod_mort) &&
            length(mod_mort) > 1 && length(mod_mort) != length(matlist)) {
          stop("mod_mort vector must be length 1 or equal to the number of demographic matrices.")
        } 
        
        if (is.list(mod_mort)) {
          if (length(mod_mort) != length(matlist)) {
            stop("mod_mort list length must match number of demographic matrices.")
          } 
          
          for (j in seq_along(mod_mort)) {
            if (!all(dim(mod_mort[[j]]) == dim(matlist[[j]]) - c(1, 0))) {
              stop(paste("mod_mort matrix at index", j, "does not match dimensions of the demographic matrix."))
            }
          }
        }
        for (i in seq_along(matlist)) {
          B <- matlist[i]
          M <- -log(B[[1]][-1, ])  # Transform survival probabilities to mortality rates
          if (is.vector(mod_mort) && length(mod_mort) == 1) {
            M <- M + mod_mort
          } else if (!is.list(mod_mort) &&
                     length(mod_mort) == length(matlist)) {
            M <- M + mod_mort[i]
          } else if (is.list(mod_mort)) {
            M <- M + mod_mort[[i]]
          }
          B[[1]][-1, ] <- exp(-M)
          matlist[i] <- B
        }
        BB <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 group_by = group_by, 
                                 lh_order = A_lh_order)
      }
## Movement modification
      if (all(!is.na(mod_move)) && is.data.frame(mod_move)) {
        matlist <- unblk.diag(MM, n_patches)
        for (i in seq_along(matlist)) {
          M <- matlist[i]
          mm <- mod_move[mod_move$stage == i, ]
          for (j in unique(mm$triangle)) {
            if (j == "lower") {
              if ("at" %in% colnames(mm)) {
                for (k in unique(mm$at)) {
                  M[lower.tri(M) & (col(M) == k)] <- 
                    M[lower.tri(M) & (col(M) == k)] * mm[mm$at == k, "d_perc"]
                  M[k, k] <- M[k, k] + (1 - sum(M[, k]))
                }
              } else if (all(c("from", "to") %in% colnames(mm))) {
                stop("Funcitnality for from - to specification coming soon.")
              }
            } else if (j == "upper") {
              if ("at" %in% colnames(mm)) {
                for (k in unique(mm$at)) {
                  M[upper.tri(M) & (col(M) == k)] <- 
                    M[upper.tri(M) & (col(M) == k)] * mm[mm$at == k, "d_perc"]
                  M[k, k] <- M[k, k] + (1 - sum(M[, k]))
                }
              } else if (all(c("from", "to") %in% colnames(mm))) {
                stop("Funcitnality for from - to specification coming soon.")
              }
            }
          }
          matlist[i] <- M
        }
        MM <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 group_by = group_by,
                                 lh_order = A_lh_order)
      }
      # ## Density-dependence
      #       for (t in 2:n_timesteps) {
      #         if (any(!is.null(ddf))){
      #           matlist <- unblk.diag(BB, n_stages)
      #           for (i in seq_along(matlist)) {
      #             B <- matlist[i]
      #             if (ddf$type == "Ricker") {
      #               B[[1]][1, ] <- dd.rec.Ricker(mat[, t - 1], ddf$a[i], ddf$b[i], theta)
      #             } 
      #             if (ddf$type == "Beverton-Holt") {
      #               B[[1]][1, ] <- dd.rec.BevertonHolt(mat[, t - 1], ddf$a[i], ddf$b[i], theta)
      #             }
      #             if (ddf$type == "logistic") {
      #               B[[1]][1, ] <- dd.growth.logistic(N = mat[c(i * n_stages - 1):c(i * n_stages), t - 1], 
      #                                                 r = ddf$r[i], B = B[[1]][1, ], K = ddf$K[i])
      #             }
      #             if (ddf$type == "ddExponential") {
      #               B[[1]][1, ] <- dd.growth.exponential(mat[, t - 1], ddf$r[i], ddf$K[i])
      #             }
      #             if (ddf$type == "general") {
      #               B[[1]][1, ] <- dd.growth.general(mat[, t - 1], ddf$r[i], ddf$K[i], ddf$theta)
      #             }
      #             matlist[i] <- B
      #           }
      #           BB <- blk.diag(matlist)
      #           A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
      #                                    group_by = group_by, lh_order = A_lh_order)
      #         }
      # ## Projection
      #         if (all(mat[, t - 1]%%1==0)) {
      #           mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
      #         } else {
      #           mat[, t] <- as.vector(A %*% mat[, t - 1]) 
      #         }
      #       }
      #       if (!is.null(rownames(n))) {
      #         rownames(mat) <- rownames(n)
      #       }
      #       colnames(mat) <- paste(1:n_timesteps)
      
      for (t in 2:n_timesteps) {
        if (!is.null(ddf)) {
          matlist <- unblk.diag(BB, n_stages)
          for (i in seq_along(matlist)) {
            B <- matlist[i]
            idx <- ((i - 1) * n_stages + 1):(i * n_stages)
            Ni <- mat[idx, t - 1]
            B[[1]][1, ] <- dd.growth.logistic(N = Ni, B = B[[1]][1, ], K = ddf$K[i],
                                              beta = ddf$beta, theta = ddf$theta)
            matlist[i] <- B
          }
          BB <- blk.diag(matlist)
          A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                   group_by = group_by, lh_order = A_lh_order)
        }
        if (anyNA(mat[, t - 1])) {
          warning(paste("NA detected at timestep", t - 1))
          print(which(is.na(mat[, t - 1])))
        }
        
        if (all(mat[, t - 1] %% 1 == 0)) {
          mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
        } else {
          mat[, t] <- as.vector(A %*% mat[, t - 1])
        }
      }
      if (!is.null(rownames(n))) {
        rownames(mat) <- rownames(n)
      }
      colnames(mat) <- paste(1:n_timesteps)
      
# LH ORDER: stages - demo -------------------------------------------------
    } else if (lh_order == "stages demo") {
      mat <- matrix(nrow = n_patches * n_stages, ncol = n_timesteps)
      tryCatch({
        mat[, 1] <- as.vector(n)
      }, error = function(e) {
        stop("Length of n and n_stages × n_patches are not equal.")
      })
## Mortality modification      
      if (any(!is.na(mod_mort))) {
        matlist <- unblk.diag(BB, n_stages)
        if (is.vector(mod_mort) &&
            length(mod_mort) > 1 && length(mod_mort) != length(matlist)) {
          stop("mod_mort vector must be length 1 or equal to the number of demographic matrices.")
        } 
        
        if (is.list(mod_mort)) {
          if (length(mod_mort) != length(matlist)) {
            stop("mod_mort list length must match number of demographic matrices.")
          } 
          
          for (j in seq_along(mod_mort)) {
            if (!all(dim(mod_mort[[j]]) == dim(matlist[[j]]) - c(1, 0))) {
              stop(paste("mod_mort matrix at index", j, "does not match dimensions of the demographic matrix."))
            }
          }
        }
        for (i in seq_along(matlist)) {
          B <- matlist[i]
          M <- -log(B[[1]][-1, ])  # Transform survival probabilities to mortality rates
          if (is.vector(mod_mort) && length(mod_mort) == 1) {
            M <- M + mod_mort
          } else if (!is.list(mod_mort) &&
                     length(mod_mort) == length(matlist)) {
            M <- M + mod_mort[i]
          } else if (is.list(mod_mort)) {
            M <- M + mod_mort[[i]]
          }
          B[[1]][-1, ] <- exp(-M)
          matlist[i] <- B
        }
        BB <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 group_by = group_by, 
                                 lh_order = A_lh_order)
      }
## Movement modification
      if (all(!is.na(mod_move)) && is.data.frame(mod_move)) {
        matlist <- unblk.diag(MM, n_patches)
        for (i in seq_along(matlist)) {
          M <- matlist[i]
          mm <- mod_move[mod_move$stage == i, ]
          for (j in unique(mm$triangle)) {
            if (j == "lower") {
              if ("at" %in% colnames(mm)) {
                for (k in unique(mm$at)) {
                  M[lower.tri(M) & (col(M) == k)] <- 
                    M[lower.tri(M) & (col(M) == k)] * mm[mm$at == k, "d_perc"]
                  M[k, k] <- M[k, k] + (1 - sum(M[, k]))
                }
              } else if (all(c("from", "to") %in% colnames(mm))) {
                stop("Funcitnality for from - to specification coming soon.")
              }
            } else if (j == "upper") {
              if ("at" %in% colnames(mm)) {
                for (k in unique(mm$at)) {
                  M[upper.tri(M) & (col(M) == k)] <- 
                    M[upper.tri(M) & (col(M) == k)] * mm[mm$at == k, "d_perc"]
                  M[k, k] <- M[k, k] + (1 - sum(M[, k]))
                }
              } else if (all(c("from", "to") %in% colnames(mm))) {
                stop("Funcitnality for from - to specification coming soon.")
              }
            }
          }
          matlist[i] <- M
        }
        MM <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 group_by = group_by,
                                 lh_order = A_lh_order)
      }
      # ## Density-dependence
      #       for (t in 2:n_timesteps) {
      #         if (any(!is.null(ddf))){
      #           matlist <- unblk.diag(BB, n_stages)
      #           for (i in seq_along(matlist)) {
      #             B <- matlist[i]
      #             if (ddf$type == "Ricker") {
      #               B[[1]][1, ] <- dd.rec.Ricker(mat[, t - 1], ddf$a[i], ddf$b[i], theta)
      #             } 
      #             if (ddf$type == "Beverton-Holt") {
      #               B[[1]][1, ] <- dd.rec.BevertonHolt(mat[, t - 1], ddf$a[i], ddf$b[i], theta)
      #             }
      #             if (ddf$type == "logistic") {
      #               B[[1]][1, ] <- dd.growth.logistic(N = mat[c(i * n_stages - 1):c(i * n_stages), t - 1], 
      #                                                 r = ddf$r[i], B = B[[1]][1, ], K = ddf$K[i])
      #             }
      #             if (ddf$type == "ddExponential") {
      #               B[[1]][1, ] <- dd.growth.exponential(mat[, t - 1], ddf$r[i], ddf$K[i])
      #             }
      #             if (ddf$type == "general") {
      #               B[[1]][1, ] <- dd.growth.general(mat[, t - 1], ddf$r[i], ddf$K[i], ddf$theta)
      #             }
      #             matlist[i] <- B
      #           }
      #           BB <- blk.diag(matlist)
      #           A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
      #                                    group_by = group_by, lh_order = A_lh_order)
      #         }
      # ## Projection
      #         if (all(mat[, t - 1]%%1==0)) {
      #           mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
      #         } else {
      #           mat[, t] <- as.vector(A %*% mat[, t - 1]) 
      #         }
      #       }
      #       if (!is.null(rownames(n))) {
      #         rownames(mat) <- rownames(n)
      #       }
      #       colnames(mat) <- paste(1:n_timesteps)
      
      for (t in 2:n_timesteps) {
        if (!is.null(ddf)) {
          matlist <- unblk.diag(BB, n_stages)
          for (i in seq_along(matlist)) {
            B <- matlist[i]
            idx <- ((i - 1) * n_stages + 1):(i * n_stages)
            Ni <- mat[idx, t - 1]
            B[[1]][1, ] <- dd.growth.logistic(N = Ni, B = B[[1]][1, ], K = ddf$K[i],
                                              beta = ddf$beta, theta = ddf$theta)
            matlist[i] <- B
          }
          BB <- blk.diag(matlist)
          A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                   group_by = group_by, lh_order = A_lh_order)
        }
        if (anyNA(mat[, t - 1])) {
          warning(paste("NA detected at timestep", t - 1))
          print(which(is.na(mat[, t - 1])))
        }
        
        if (all(mat[, t - 1] %% 1 == 0)) {
          mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
        } else {
          mat[, t] <- as.vector(A %*% mat[, t - 1])
        }
      }
      if (!is.null(rownames(n))) {
        rownames(mat) <- rownames(n)
      }
      colnames(mat) <- paste(1:n_timesteps)

# LH ORDER: stages - move -------------------------------------------------
    } else if (lh_order == "stages move") {
      mat <- matrix(nrow = n_patches * n_stages, ncol = n_timesteps)
      tryCatch({
        mat[, 1] <- as.vector(n)
      }, error = function(e) {
        stop("Length of n and n_stages × n_patches are not equal.")
      })
## Mortality modificaiton      
      if (any(!is.na(mod_mort))) {
        matlist <- unblk.diag(BB, n_stages)
        if (is.vector(mod_mort) &&
            length(mod_mort) > 1 && length(mod_mort) != length(matlist)) {
          stop("mod_mort vector must be length 1 or equal to the number of demographic matrices.")
        } 
        
        if (is.list(mod_mort)) {
          if (length(mod_mort) != length(matlist)) {
            stop("mod_mort list length must match number of demographic matrices.")
          } 
          
          for (j in seq_along(mod_mort)) {
            if (!all(dim(mod_mort[[j]]) == dim(matlist[[j]]) - c(1, 0))) {
              stop(paste("mod_mort matrix at index", j, "does not match dimensions of the demographic matrix."))
            }
          }
        }
        for (i in seq_along(matlist)) {
          B <- matlist[i]
          M <- -log(B[[1]][-1, ])  # Transform survival probabilities to mortality rates
          if (is.vector(mod_mort) && length(mod_mort) == 1) {
            M <- M + mod_mort
          } else if (!is.list(mod_mort) &&
                     length(mod_mort) == length(matlist)) {
            M <- M + mod_mort[i]
          } else if (is.list(mod_mort)) {
            M <- M + mod_mort[[i]]
          }
          B[[1]][-1, ] <- exp(-M)
          matlist[i] <- B
        }
        BB <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 group_by = group_by, 
                                 lh_order = A_lh_order)
      }
## Movement modification
      if (all(!is.na(mod_move)) && is.data.frame(mod_move)) {
        matlist <- unblk.diag(MM, n_patches)
        for (i in seq_along(matlist)) {
          M <- matlist[i]
          mm <- mod_move[mod_move$stage == i, ]
          for (j in unique(mm$triangle)) {
            if (j == "lower") {
              if ("at" %in% colnames(mm)) {
                for (k in unique(mm$at)) {
                  M[lower.tri(M) & (col(M) == k)] <- 
                    M[lower.tri(M) & (col(M) == k)] * mm[mm$at == k, "d_perc"]
                  M[k, k] <- M[k, k] + (1 - sum(M[, k]))
                }
              } else if (all(c("from", "to") %in% colnames(mm))) {
                stop("Funcitnality for from - to specification coming soon.")
              }
            } else if (j == "upper") {
              if ("at" %in% colnames(mm)) {
                for (k in unique(mm$at)) {
                  M[upper.tri(M) & (col(M) == k)] <- 
                    M[upper.tri(M) & (col(M) == k)] * mm[mm$at == k, "d_perc"]
                  M[k, k] <- M[k, k] + (1 - sum(M[, k]))
                }
              } else if (all(c("from", "to") %in% colnames(mm))) {
                stop("Funcitnality for from - to specification coming soon.")
              }
            }
          }
          matlist[i] <- M
        }
        MM <- blk.diag(matlist)
        A <- spmm.project.matrix(P = P, BB = BB, MM = MM,
                                 group_by = group_by,
                                 lh_order = A_lh_order)
      }
      # ## Density-dependence
      #       for (t in 2:n_timesteps) {
      #         if (any(!is.null(ddf))){
      #           matlist <- unblk.diag(BB, n_stages)
      #           for (i in seq_along(matlist)) {
      #             B <- matlist[i]
      #             if (ddf$type == "Ricker") {
      #               B[[1]][1, ] <- dd.rec.Ricker(mat[, t - 1], ddf$a[i], ddf$b[i], theta)
      #             } 
      #             if (ddf$type == "Beverton-Holt") {
      #               B[[1]][1, ] <- dd.rec.BevertonHolt(mat[, t - 1], ddf$a[i], ddf$b[i], theta)
      #             }
      #             if (ddf$type == "logistic") {
      #               B[[1]][1, ] <- dd.growth.logistic(N = mat[c(i * n_stages - 1):c(i * n_stages), t - 1], 
      #                                                 r = ddf$r[i], B = B[[1]][1, ], K = ddf$K[i])
      #             }
      #             if (ddf$type == "ddExponential") {
      #               B[[1]][1, ] <- dd.growth.exponential(mat[, t - 1], ddf$r[i], ddf$K[i])
      #             }
      #             if (ddf$type == "general") {
      #               B[[1]][1, ] <- dd.growth.general(mat[, t - 1], ddf$r[i], ddf$K[i], ddf$theta)
      #             }
      #             matlist[i] <- B
      #           }
      #           BB <- blk.diag(matlist)
      #           A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
      #                                    group_by = group_by, lh_order = A_lh_order)
      #         }
      # ## Projection
      #         if (all(mat[, t - 1]%%1==0)) {
      #           mat[, t] <- floor(as.vector(A %*% mat[, t - 1]))
      #         } else {
      #           mat[, t] <- as.vector(A %*% mat[, t - 1]) 
      #         }
      #       }
      #       if (!is.null(rownames(n))) {
      #         rownames(mat) <- rownames(n)
      #       }
      #       colnames(mat) <- paste(1:n_timesteps)
      
      for (t in 2:n_timesteps) {
        if (!is.null(ddf)) {
          matlist <- unblk.diag(BB, n_stages)
          for (i in seq_along(matlist)) {
            B <- matlist[i]
            idx <- ((i - 1) * n_stages + 1):(i * n_stages)
            Ni <- mat[idx, t - 1]
            B[[1]][1, ] <- dd.growth.logistic(N = Ni, B = B[[1]][1, ], K = ddf$K[i],
                                              beta = ddf$beta, theta = ddf$theta)
            matlist[i] <- B
          }
          BB <- blk.diag(matlist)
          A <- spmm.project.matrix(P = P, BB = BB, MM = MM, 
                                   group_by = group_by, lh_order = A_lh_order)
        }
        if (anyNA(mat[, t - 1])) {
          warning(paste("NA detected at timestep", t - 1))
          print(which(is.na(mat[, t - 1])))
        }
        
        if (all(mat[, t - 1] %% 1 == 0)) {
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

# Package output ----------------------------------------------------------
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
        "A projection matrix. The arg mod_move is currently ignored; please modify manually."
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
