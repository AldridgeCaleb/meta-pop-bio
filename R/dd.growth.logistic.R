##' @title Density dependent logistic growth
##' 
##' @description
##' Adjusts per capita fecundity based on current total population size at *t*
##' relative to carrying capacity using a theta-logistic-like formulation. 
##' Allows enhancement at low density and compensation at high density.
##'
##' @param N Numeric vector. Abundance or biomass at time *t* for each stage or group.
##' @param B Numeric vector. Baseline or average per capita fecundity rates 
##' for each stage. Typically, only mature stages are nonzero.
##' @param K Numeric scalar. Carrying capacity, representing max abundance or biomass 
##' in the population, beyond which fecundity begins to decline.
##' @param beta Numeric scalar. Density feedback strength. 
##' \code{beta = 0} implies no density dependence. 
##' \code{beta > 0} enables enhanced reproduction at low density.
##' @param theta Numeric scalar. Shape parameter controlling the steepness of the 
##' response. \code{theta = 1} is standard logistic; \code{theta > 1} produces 
##' stronger decline at high density.
##'
##' @note
##' This theta-logistic formulation allows low-density enhancement (if \code{beta > 0}), 
##' mimicking undercompensatory or compensatory responses. 
##'
##' @return A numeric vector of modified per capita fecundity rates (same length as \code{N}), 
##' adjusted for density dependence.
##' 
##' @export
dd.growth.logistic <- function(N, B, K, beta = 0, theta = 1) {
  if (is.null(K) || K <= 0) stop("K must be provided and > 0.")
  if (is.null(B)) stop("B must be provided.")
  
  mat_idx <- which(B > 0)
  N_total <- sum(N)
  
  # Compute ratio and protect against non-finite values
  ratio <- N_total / K
  if (!is.finite(ratio)) ratio <- 1e6  # large fallback value
  
  # Density-dependent enhancement possible when beta > 0
  raw_modifier <- 1 + beta * (1 - ratio^theta)
  if (!is.finite(raw_modifier) || is.nan(raw_modifier)) raw_modifier <- 0
  density_modifier <- max(density_modifier, 0)  # truncate negative values
  
  # Apply modifier to baseline per capita fertility
  out <- rep(0, length(N))
  out[mat_idx] <- B[mat_idx] * density_modifier
  
  return(out)
}


