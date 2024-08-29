##' @title Density-dependent Beverton-Holt recruitment function
##'
##' @param N Number of individuals or biomass in time t. Can (should) be a 
##' vector structure according to group_by.
##' @param a Density-independent recruitment rate.
##' @param b Rate of change with increased density.
##' @param theta Degree of depensation (Allee effect).
##'
##' @note
##' Equations expressed as recruit per capita.
##' 
##' @export
dd.rec.BevertonHolt <- function(N, a, b, theta = NA) {
  N <- sum(N)  # Sum of all individuals or biomass
  if (!is.na(theta) && theta != 1) {
    # Depensation Beverton-Holt with Allee effect
    return((a * N^(theta - 1)) / (1 + b * N^theta))
  } else {
    # Standard Beverton-Holt
    return(a / (1 + b * N))
  }
}
