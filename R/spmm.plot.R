#' @title plots of spatial matrix population model projections
#'
#' @description
#' Produces plots of deterministic future population sizes of stages s in
#' patches p. Depending on structure of `n`, plot frames represent patches or
#' stages and lines the opposite.
#'
#' @param projections A `matrix` with `nrow` n_stages × n_patches and `ncol` of
#' n_timesteps (see `spmm.project` for more detail). If prerequisite steps,
#' i.e., `vec.perm`, `blk.diag`, `spmm.project.matrix`, and `spmm.project`, have been
#' specified correctly and correspond to structure of `n` (and expectations).
#' @param ylabs Y-axis label for plots.
#' @param xlabs X-axis label for plots.
#' @param stage_names Names of stages, ages, classes, etc.
#' @param patch_names Names of patches, units, pools, etc.
#' @param ylim_max Choice of either "overall" or "patch" that controls the ylim
#' for each plot.
#'
#' @note
#' As with `spmm.project` ensure that the structural type of population vector
#' `n` and projection matrix `A` are the same. Otherwise, projections may produce
#' incorrect values!
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
#' # Peregrine falcon example from Hunter and Caswell (2005), Ecological Modelling 
#' # 188(2005):15--21. Data from Wootton and Bell (1992). Continues example from 
#' # `spmm.project`.
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
#' # Initial stages within patches (patch group_by)
#' n <- c(
#'   50, 22,  # Northern patch adults and juveniles
#'   40, 17   # Southern patch adults and juveniles
#' )
#' comment(n) <- "patches"  # vec comment attr for group_by
#'
#' # Number of time steps to project into the future
#' n_timesteps <- 50
#'
#' # Project spatial matrix population model
#' projs <- spmm.project(n, A, n_timesteps, n_stages, n_patches)
#'
#' # Plot projections
#' spmm.plot(projs)
#'
#' @export
spmm.plot <- function(projections, ylabs = NA, xlabs = NA, 
                          stage_names = NA, patch_names = NA,
                      ylim_max = "overall") {
  comments <- comment(projections)
  group_by <- strsplit(comments, " +")[[1]][1]
  if (group_by == "patches") {
    n_stages <- as.numeric(strsplit(comments, " +")[[2]][1])
    n_patches <- as.numeric(strsplit(comments, " +")[[3]][1])
    
    graphics::par(
      mfrow = c(min(round(n_patches / 2, 0), 3), 2),
      mar = c(5, 5, 1.5, 0.5),
      oma = rep(0.5, 4)
    )
    
    starts <- seq(1, dim(projections)[1], by = n_stages)
    ends <- c(starts - 1, dim(projections)[1])[-1]
    # throw error if starts and ends != lengths()
    for (i in seq_along(starts)) {
      idx <- starts[i]:ends[i]
      if (ylim.max == "patch") {
        ylim_vals <- c(0, round(max(projections[idx, , drop = FALSE]) + 1, -1))
      } else {
        ylim_vals <- c(0, round(max(projections) + 1, -1))
      }
      graphics::matplot(
        t(projections)[, idx],
        type = 'b',
        pch = 16,
        # col = c("black", "black"),
        ylim = ylim_vals,
        ylab = ylabs,
        xlab = xlabs,
        main = paste("Patch :", patch_names[i])
      )
    }
    if (!is.na(stage_names)[1]) {
      graphics::plot.new()
      graphics::legend("center",
             stage_names,
             pch = 16,
             col = 1:length(stage_names))
    }
    
  } else if (group_by == "stages") {
    n_patches <- as.numeric(strsplit(comments, " +")[[2]][1])
    n_stages <- as.numeric(strsplit(comments, " +")[[3]][1])
    
    graphics::par(
      mfrow = c(min(round(n_patches / 2, 0), 3), 2),
      mar = c(5, 5, 1.5, 0.5),
      oma = rep(0.5, 4)
    )
    
    starts <- seq(1, dim(projections)[1], by = n_patches)
    ends <- c(starts - 1, dim(projections)[1])[-1]
    # throw error if starts and ends != lengths()
    for (i in seq_along(starts)) {
      idx <- starts[i]:ends[i]
      if (ylim.max == "stage") {
        ylim_vals <- c(0, round(max(projections[idx, , drop = FALSE]) + 1, -1))
      } else {
        ylim_vals <- c(0, round(max(projections) + 1, -1))
      }
      graphics::matplot(
        t(projections)[, idx],
        type = 'b',
        pch = 16,
        # col = c("black", "black"),
        ylim = ylim_vals,
        ylab = ylabs,
        xlab = xlabs,
        main = paste("Stage :", stage_names[i])
      )
    }
    if (!is.na(patch_names)[1]) {
      graphics::plot.new()
      graphics::legend("center",
             patch_names,
             pch = 16,
             col = 1:length(patch_names))
    }
  }
}
