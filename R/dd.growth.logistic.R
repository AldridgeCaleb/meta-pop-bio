##' @title Density dependent logistic growth
##' 
##' @description
##' Modifies growth of population in time *t + 1* by the population size in 
##' time *t*.
##'
##' @param N Number of individuals or biomass in time t. Can (should) be a 
##' vector structure according to group_by.
##' @param r Intrinsic rate of increase or baseline growth rate.
##' @param K Carrying capacity. Can be a vector for each stage and can vary by
##' patches. 
##'
##' @note
##' This produces similar results as `dd.growth.exponential`.
##' 
##' @export
dd.growth.logistic <- function(N, r = NULL, B = NULL, K) {
  if (is.null(r)) {
    if (is.null(B) || all(is.na(B))) {
      stop("Both 'r' and valid 'B' must be provided for logistic density dependence.")
    }
    r <- N * B
  }
  N_sum <- sum(N)
  return(r * (1 - N_sum/K))
}
