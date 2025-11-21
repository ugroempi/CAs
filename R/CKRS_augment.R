#' Function to augment a CA with Eqn (5) of CKRS
#'
#' increase the number of levels of a CA by 1 in a greedy fashion
#'
#' @rdname CKRS_augment
#'
#' @aliases CKRS_augmentCA
#' @aliases CKRS_augment
#' @aliases N_CKRS_augmentCA
#'
#' @usage CKRS_augmentCA(t, k, v, ...)
#' @usage CKRS_augment(D, start0=TRUE, ...)
#' @usage N_CKRS_augmentCA(t, k, v, ...)
#'
#' @param t integer: the target strength
#' @param k integer: the target number of columns
#' @param v integer: the number of levels
#' @param D a strength t CA with k columns in v-1 levels, coded as integers
#'           0, ..., v-2 (or 1, ..., v-1)
#' @param start0 logical: Do integer values start with 0 ?
#' @param ... further arguments to function \code{\link{bestN}}
#'       (function \code{N_CKRS_augmentCA})
#'
#' @details
#' \code{CKRS_augmentCA} is the interface function that uses the workhorse
#' function \code{CKRS_augment}.
#'
#' @returns Function \code{CKRS_augmentCA} returns a matrix of class \code{ca}
#' with attributes, which is a strength t covering array.\cr
#' Function \code{CKRS_augment} returns a matrix without class and attributes,\cr
#' function \code{N_CKRS_augmentCA} returns the run size of the CA.
#'
#' Note that the run size may depend on the availability of an internet
#' connection for cases for which the best strength t CA for \code{k}
#' columns with v-1 levels is
#' in the Dwyer data base or the NIST library only
#' (these might exist - has not been checked).
#'
#' @references CKRS (2010)
#'
#' @examples
#' dim(CKRS_augment(SCA_Busht(5,3))) ## not competitive
#' bestN(3,6,6)
#'
#' @export
CKRS_augment <- function(D, start0=TRUE, ...){
  if (!is.matrix(D)) D <- as.matrix(D)
  ll <- levels.no(D)
  stopifnot(length(unique(ll))==1)
  v <- ll[1]
  k <- ncol(D)
  if (start0) stopifnot(min(D)==0) else stopifnot(min(D)==1)
  if (!start0) D <- D - 1
  for (i in 1:k){
    tab <- table(D[,i])
    pick <- which.min(tab)
    lev <- as.numeric(names(tab)[pick])
    addo <- D[which(D[,i]==lev),,drop=FALSE]
    addo[,i] <- v
    D <- rbind(D, addo)
  }

  if (!start0) D <- D + 1
  return(D)
}

#' @export
CKRS_augmentCA <- function(t, k, v, ...){
  if (v==2) stop("augmentation does not work for v=2")
  Call <- sys.call()
  D <- bestCA(t, k, v-1)
  aus <- CKRS_augment(D)
  class(aus) <- c("ca", class(aus))
  attr(aus, "origin") <- "augmentation according to Eqn (5), CKRS"
  attr(aus, "Call") <- Call
  attr(aus, "t") <- t
  aus
}

#' @export
N_CKRS_augmentCA <- function(t, k, v, ...){
  if (v==2) return(NA)
  M <- bestN(t, k, v-1)
  for (i in 1:k)
    M <- floor(M*v/(v-1))
  unname(M)
}
