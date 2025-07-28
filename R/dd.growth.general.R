##' @title density dependent growth generalized
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
##' @param theta Shape parameter that modifies the growth form from *t* to *t+1*.
##' If *theta* = 1, then it takes a logitic shape (Beverton-Holt-like); if > 1, 
##' it becomes more concave; if < 1, it becomes more convex; and if *theta* = 0, 
##' then it takes an density-dependent exponential shape (Ricker-like). 
##' 
##' @note
##' This produces very similar results as `dd.growth.logistic`.
##'
##' @export
dd.growth.general <- function(N, r, K, theta = 1) {
  N <- sum(N)
  return((1 + r * (1 - N / K)) ^ theta)
}
