##' @title ensity dependent logistic growth
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
##' This produces very similar results as `dd.growth.exponential`.
##' 
##' @export
dd.growth.logistic <- function(N, r, K) {
  N <- sum(N)
  return(1 + (r * (1 - N/K)))
}
