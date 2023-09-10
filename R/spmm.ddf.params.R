#' @title define density dependence function parameters
#' 
#' @param f_type Either Ricker or Beverton-Holt type density dependence on 
#' recruitment.
#' @param s_type Currently implments only logistic type density dependence on 
#' survival.
#' @param stages Number of stages being modeled.
#' @param P Vec-permutation matrix.
#' @param BB Block diagonal demographic matrix.
#' @param MM Block diagonal movement (dispersal) matrix.
#' @param r Baseline recruitment.
#' @param K Carrying capacity.
#' 
#' @export
spmm.ddf.params <- function(f_type = NA, s_type = NA, stages = NA,
                     P = NA, BB = NA, MM = NA, 
                     r = NA, K = NA) {
  return(list(f_type, s_type, stages,
              P, BB, MM, r, K))
} 
