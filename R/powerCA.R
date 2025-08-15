#' Power constructions
#'
#' Functions to construct a strength t CA in with k columns at v levels
#' based on the power construction by Colbourn and Torres-Jimenez (2010)
#'
#' @aliases powerCA
#' @aliases N_powerCT
#' @aliases k_powerCT
#'
#' @usage powerCA(t, k, v, type="CT", ...)
#' @usage N_powerCT(t, k, v, type="CT", ...)
#' @usage k_powerCT(t, N, v, type="CT", ...)
#'
#' @param t integer, strength of the CA
#' @param k integer, number of columns of the CA
#' @param v integer, number of levels of the CA
#' @param N integer, number of runs of the CA
#' @param type character string, currently only \code{CT} is implemented
#' @param ... currently not used
#'
#' @section Details:
#' The function \code{powerCA} implements constructions according to
#' Colbourn and Torres-Jimenez (2010); it is based on the data.frame object
#' \code{\link{powerCTcat}} that contains scenarii for which the power construction
#' with ideal ingredients yields the smallest possible CAs according to
#' \code{\link{colbournBigFrame}} (and for which the construction string in
#' \code{\link{colbournBigFrame}} is understood).
#'
#' The actual arrays that can be generated with this package may be larger than
#' the best possible ones if the perfect ingredients are not (yet) available;
#' they may also be somewhat smaller, because this package makes use of the
#' maximum number of constant rows of the ingredient CAs where possible, whereas
#' the numbers in \code{\link{colbournBigFrame}} do not always do that (they do it where
#' the Source entry contains "c" next to the respective ingredient).
#'
#' The construction uses functions \code{\link{Tred}}, \code{\link{createDHF}}, and
#' \code{\link{DHHF2CA}}, as well as CA construction functions as indicated in
#' \code{\link{powerCTcat}}.
#'
#' @returns Function \code{powerCA}
#' returns a strength \code{t} CA in \code{N_powerCA(t, k, v)} runs
#' and \code{k} columns. \code{k_powerCA(t, N, v)} returns the maximum number of columns
#' that can be achieved with \code{N} runs when requesting strength \code{t} at \code{v} levels.
#'
#' @references Colbourn and Torres-Jimenez (2010)
#'
#' @examples
#' Ns(3, 25, 2)  ## not a good case, but small enough for checking coverage
#' D <- powerCA(3, 25, 2)
#' dim(D)
#' coverage(D, 3)
#' Ns(3, 25, 2) ## too large for so few columns
#'
#' ## a larger example
#' Ns(6, 25, 2) ## way too large for so few columns
#' ks(6, 5321, 2) ## would cover many more columns
#' dim(powerCA(6,5000,2))
#' Ns(6, 5000, 2) ## eCAN is PowerCZ
#'

#' @export
powerCA <- function(t, k, v, type="CT", ...){
  Call <- sys.call()
  pick <- which(powerCTcat$t >= t &
                  powerCTcat$k >= k &
                  powerCTcat$v == v )
  if (length(pick)==0) stop("no construction for this setting")
  hilf <- powerCTcat[pick,,drop=FALSE]
  r <- which.min(hilf$N)
  ## provide new t, which can be larger than requested t
  t <- hilf$t[r]
  oaDHF <- eval(parse(text=hilf$OAforDHF[r]))
  toa <- hilf$expon[r]
  w1 <- hilf$w1[r]
  DHF <- createDHF(oaDHF, toa, t, v)
  M <- nrow(DHF)
  stopifnot(M==hilf$M[r])

  ## first OA
  N1 <- hilf$N1[r]
  D1 <- eval(parse(text=labelToCode(hilf$constr1[r], t, w1, v)))[,1:w1]
  if (N1 > 250) D1 <- maxconstant(D1, one_is_enough = TRUE) else D1 <- maxconstant(D1)

  OAlist <- list(D1)

  u2 <- hilf$u2[r]; u3 <- hilf$u3[r]; u4 <- hilf$u4[r]

  ## handle heterogeneity
  if (u2 > 0){
    w2 <- hilf$w2[r]
    N2 <- hilf$N2[r]
    D2=eval(parse(text=labelToCode(hilf$constr2[r], t, w2, v)))[,1:w2]
    if (N2 > 250) D2 <- maxconstant(D2, one_is_enough = TRUE) else D2 <- maxconstant(D2)
    OAlist <- append(OAlist, list(D2))
    if (u3 > 0){
      w3 <- hilf$w3[r]
      N3 <- hilf$N3[r]
      D3 <- eval(parse(text=labelToCode(hilf$constr3[r], t, w3, v)))[,1:w3]
      if (N3 > 250) D3 <- maxconstant(D3, one_is_enough = TRUE) else D3 <- maxconstant(D3)
      OAlist <- append(OAlist, list(D3))
      if (u4 > 0){
        w4 <- hilf$w4[r]
        N4 <- hilf$N4[r]
        D4 <- eval(parse(text=labelToCode(hilf$constr4[r], t, w4, v)))[,1:w4]
        if (N4 > 250) D4 <- maxconstant(D4, one_is_enough = TRUE) else D4 <- maxconstant(D4)
        OAlist <- append(OAlist, list(D4))
        rs <- (M-u4+1):M
        red <- w1 - w4
        DHF <- Tred(DHF, rs, red)
      }
      rs <- (M-u4-u3+1):(M-u4)
      red <- w1-w3
      DHF <- Tred(DHF, rs, red)
    }
    rs <- (M-u4-u3-u2+1):(M-u4-u3)
    red <- w1-w2
    DHF <- Tred(DHF, rs, red)
  }
  aus <- DHHF2CA(DHF, OAlist, v=v, one_is_enough=TRUE)[, 1:k]
  class(aus) <- c("ca", class(aus))
  attr(aus, "Call") <- Call
  attr(aus, "t") <- t
  attr(aus, "origin") <- "Power construction of Colbourn and Torres-Jimenez 2010"
  attr(aus, "date") <- Sys.Date()
  attr(aus, "CAs-version") <- packageVersion(pkg="CAs")
  if (any(is.na(aus))) attr(aus, "flexible") <- list(value=NA, profile=flexprofile(aus))
  aus
}

#' @export
N_powerCT <- function(t,k,v, type="CT", ...){
  hilf <- powerCTcat[powerCTcat[,"t"]>=t & powerCTcat[,"k"]>=k & powerCTcat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(min(hilf[,"N"]))
}

#' @export
k_powerCT <- function(t,N,v, type="CT", ...){
  hilf <- powerCTcat[powerCTcat[,"t"]>=t & powerCTcat[,"N"]<=N & powerCTcat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(max(hilf[,"k"]))
}

#
# ## is this worth implementing?
# squareCA <- function(D, t, check=FALSE, ...){
#   ## Thm. 3.2 of Chateauneuf, Colbourn and Kreher
#   stopifnot(is.matrix(D))
#   k <- ncol(D)
#   stopifnot(is.numeric(t))
#   stopifnot(t>=2)
#   stopifnot(t%%1==0)
#   if (!primes::gcd(factorial(choose(t,2)), k)==1) stop("this construction only works, if factorial(choose(t,2)) and k have no common primes")
#   stopifnot(is.numeric(D))
#   stopifnot(min(D)%%1==0)
#   stopifnot(all(D>=0))
#   start0 <- TRUE
#   if (min(D=1)){
#     start0 <- FALSE
#     D <- D - 1
#   }
#   if (check){
#     stopifnot(all(coverage(D,t)==1))
#   }
#   ## need a difference matrix D(k, choose(t,2)+1, 1)
#
# }
