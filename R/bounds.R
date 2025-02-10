#' functions for bounds for CA sizes
#'
#' The function(s) provide bounds for CA sizes.
#'
#' @rdname bounds
#'
#' @aliases lower_bound_CA
#'
#' lower_bound_CA(t, k, v)
#'
#' @param t requested strength of CA (positive integer number >= 2)
#' @param k number of columns, must be at least \code{t}
#' @param v number of levels for each column:
#' a single integer number or a vector of integer numbers of length \code{k}
#'
#' @returns \code{lower_bound_CA} returns a positive number
#'
#' @section Details:
#' Function \code{lower_bound_CA} provides the trivial lower bound which is the
#' largest product of \code{t} numbers of levels.
#'
#' @examples
#' lower_bound_CA(3, 75, 2)  ## 2^3
#' lower_bound_CA(3, 75, c(rep(2,72), 3, 5, 8)) ## 3*5*8
#'

#' @export
## lower bound based on size of minimum t-tuple
lower_bound_CA <- function(t, k, v){
    stopifnot(k >= t)
    if (length(v)==1) return(v^t) else{
      stopifnot(length(v)==k)
      stopifnot(all(v%%1==0))
      return(prod(sort(v, decreasing=TRUE)[1:t]))
    }
}

## obtain lower bound for balanced strength 2 CA
## not exported because not perceived as useful
lower_t2_balanced <- function(k,v){
  ## based on Kokkala et al. arXiv preprint for 2020 paper
  ##    "On the Structure of Small Strength-2 Covering Arrays"
  ##    for uniform strength 2 CAs all of whose individual columns
  ##    are maximally balanced w.r.t. their level distributions
  ## this does not seem very convincing, doesn't change with k!
    a <- 2*k^2 - 6*k + 4*v
    b <- (2*v - 3)*k^2 + (-2*v + 5)*k + (-4*v + 2)
    D <- k^4 + (8*v^2 - 16*v + 2)*k^3 +
      (-8*v^3 + 24*v - 3)*k^2 + (8*v^3 - 8*v^2 - 8*v - 4)*k + 4
    max(v^2, ceiling(v*(b+sqrt(D))/a))
}
