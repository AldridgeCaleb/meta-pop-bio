##' @name dd.funs
##' @aliases dd.surv.logistic
##' @aliases dd.rec.Ricker
##' @aliases dd.rec.BevertonHolt
##'
##' @title density-dependent recruitment functions
##'
##' @param N Number of individuals or biomass in time t. Can (should) be a 
##' vector structure according to group_by.
##' @param r Baseline or maximum recruitment rate.
##' @param K Carrying capacity. Can be a vector for each stage and can vary by
##' patches. 
##'
##' @note \code{dd.funs} is a generic name for the functions documented.
##' \cr
##' If called, \code{dd.funs} returns its own arguments.
##'
##' @rdname dd.funs
##' @export
dd.funs <- function(N, r, K) {identity(c(N, r, K))}
##'
##' @rdname dd.surv.logistic
##' @return \code{dd.surv.logistic(N, r, K)}
##' @examples
##' dd.surv.logistic(100, 0.62, 200)
##' @export
dd.surv.logistic <- function(N, r, K) {
  N <- sum(N)
  return(1 - (N / K))
}
##'
##' @rdname dd.rec.Ricker
##' @return \code{dd.rec.Ricker(N, r, K)}
##' @examples
##' dd.rec.Ricker(100, 0.62, 200)
##' @export
dd.rec.Ricker <- function(N, r, K) {
  N <- sum(N)
  return(exp(r * (1 - (N/K))))
}
##'
##' @rdname dd.rec.BevertonHolt
##' @return \code{dd.rec.BevertonHolt(N, r, K)}
##' @examples
##' dd.rec.BevertonHolt(100, 0.62, 200)
##' @export
dd.rec.BevertonHolt <- function(N, r, K) {
  N <- sum(N)
  return(r / (1 + (r - 1) * N / K))
}
