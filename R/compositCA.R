#' Obtain a CA for non-prime number of levels by composition
#'
#' Function to construct a strength t CA in k columns at v levels, v non-prime, from two strength t CAs in fewer levels
#'
#' @rdname compositCA
#'
#' @aliases compositCA
#' @aliases N_compositCA
#'
#' @usage compositCA(t, k, v, ...)
#' @usage N_compositCA(t, k, v, ...)
#'
#' @param t positive integer: requested strength
#' @param k positive integer: requested number of columns
#' @param v positive integer, non-prime: requested number of levels for each column
#' @param ... currently not used
#'
#' @section Details:
#' Function \code{compositCA} factorizes \code{v} into all variants of two factors,
#' picks the pair that produces the smallest run size, and obtains a CA in \code{v} levels
#' by composition from two CAs with those numbers of levels, using function \code{\link{crossCAs}}.
#'
#' It uses function \code{\link{bestCA}} for obtaining the ingredient CAs. The result is rarely optimal.
#'
#' @returns Function \code{compositCA} yields a matrix of class \code{ca}), with levels from 0 to \code{v-1}.
#'
#' @examples
#' ################################################
#' ## compositCA
#' ################################################
#' D <- compositCA(3, 12, 9)
#' dim(D)
#' Ns(3, 12, 9) ## not competitive
#'
#' D <- compositCA(4, 12, 4)
#' dim(D)
#' Ns(4, 12, 4) ## best implemented
#' eCAN(4, 12, 4)  ## this has not been implemented
#'
#' Ns(3, 4, 12) ## best (orthogonal array)
#'
#' Ns(4, 24, 12) ## slightly worse than fuseBushtCA
#'
#' Ns(4, 12, 8)  ## best implemented
#'
#' Ns(4, 12, 24) ## much worse than fuseBushtCA,
#'               ## which almost yields eCAN
#'

#' @export
compositCA <- function(t, k, v, ...){
  Call <- sys.call()
  stopifnot(is.numeric(t), is.numeric(k), is.numeric(v))
  hilf <- conf.design::factorize(v)
  ## primes are not treated by composition
  if (length(hilf)==1) return(bestCA(t, k, v, ...))
  ## now non-prime

  divisors <- numeric(0)
  for (i in 2:sqrt(v)){
    if (v%%i==0L) divisors <- c(divisors, i)
  }

  ## identify best pair
  N_erg <- sapply(divisors, function(obj) bestN(t, k, obj)*bestN(t, k, v%/%obj))
  factor1 <- divisors[which.min(N_erg)]

  aus <- crossCAs(bestCA(t, k, factor1), bestCA(t, k, v%/%factor1))
  attr(aus, "t") <- t
  attr(aus, "Call") <- Call
   aus
}

#' @export
N_compositCA <- function(t, k, v, ...){
  hilf <- conf.design::factorize(v)
  ## primes and prime powers are not treated by composition
  if (length(hilf)==1) return(NA)
  ## now non-prime
  divisors <- numeric(0)
  for (i in 2:sqrt(v)){
    if (v%%i==0L) divisors <- c(divisors, i)
  }
  N_erg <- sapply(divisors, function(obj) bestN(t, k, obj)*bestN(t, k, v%/%obj))
  return(unname(min(N_erg)))
}
