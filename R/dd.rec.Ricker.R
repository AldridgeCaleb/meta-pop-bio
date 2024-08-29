##' @title Density-dependent Ricker recruitment function
##'
##' @param N Number of individuals or biomass in time t. Can (should) be a 
##' vector structure according to group_by.
##' @param a Density-independent recruitment rate.
##' @param b Rate of change with increased density.
##' @param theta Degree of depensation (Allee effect). Default is NA.
##'
##'##' @note
##' Equations expressed as recruit per capita.
##' 
##' @export
dd.rec.Ricker <- function(N, a, b, theta = NA) {
  N <- sum(N)  # Sum of all individuals or biomass
  if (!is.na(theta) && theta != 1) {
    # Depensation Ricker with Allee effect
    return(a * N^(theta - 1) * exp(-b * N))
  } else {
    # Standard Ricker
    return(a * exp(-b * N))
  }
}
