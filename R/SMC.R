## code for creating a CA from the CPHF
#' Making CA from SCPHF
#'
#' according to the process outlined in Sherwood, Martirosyan and Colbourn (2006).
#' It is assumed that the SCPHF refers to all non-reduced permutation vectors
#' or to all reduced permutation vectors, but not to a mix of both.
#'
#' @rdname SMC
#'
#' @aliases smcCA
#' @aliases SMC
#' @aliases N_smcCA
#' @aliases k_smcCA
#'
#' @usage smcCA(t, k, v, start0 = TRUE, ...)
#' @usage SMC(cphf, q, t, ...)
#' @usage N_smcCA(t, k, v)
#' @usage k_smcCA(t, N, v)
#'
#' @param t integer: target strength
#' @param k integer: requested number of columns
#' @param v integer: number of levels, must be a prime power
#' @param start0 logical: should the levels of the output array start with 0 (otherwise with 1)?
#' @param cphf an SCPHF
#' @param q integer: prime power
#' @param N integer: affordable run size
#' @param ... currently not used
#'
#' @returns Function \code{smcCA} returns a matrix of class \code{ca} with attributes,
#' which is a covering array of strength 3 or 4.\cr
#' Function \code{SMC} returns a matrix with entries in 0,...,v-1, which is a covering
#' array of strength \code{t}, if a valid \code{cphf} for strength \code{t} is provided;
#' the function is not meant for direct use, but is the workhorse behind \code{smcCA}.
#'
#' Functions \code{N_smcCA} and \code{k_smcCA} return the required run size or the achievable number of columns, respectively.
#'
#' @section Details:
#' Function \code{SMC} uses a covering perfect hash family of a specific form (SCPHF)
#' that works with length \code{v^t} permutation vectors for \code{v} levels
#' and strength \code{t}. Function \code{smcCA} implements CAs from this
#' construction for the SCPHFs from Sherwood et al. (2006).
#' The construction relies on a Galois field GF(\code{v}).
#' The SCPHFs for \code{v=8} or \code{v=9} require a different Galois field
#' than otherwise used in this package (created with internal function \code{mygf}).
#'
#' The data for the construction are in \code{\link{SMC_CPHFs}} and \code{\link{SMCcat}},
#' respectively.
#'
#' @references Sherwood, Martirosyan and Colbourn (2006)
#'
#' @examples
#' # applying the construction for small SCPHF from the Sherwood et al. (2006) paper
#'
#' # strength 3
#' cphf <- rbind(0:15,
#'               c(0,1,4,5,2,3,6,7,9,8,13,12,11,10,15,14))
#' D <- SMC(cphf, q=4, t=3)
#' dim(D)
#' coverage(D,3)
#' eCAN(3,16,4)   ## size matches eCAN
#'
#' # strength 4
#' cphf <- rbind(0:7,
#'               c(0, 1, 2, 4, 3, 6, 7, 5))
#' D <- SMC(cphf, q=2, t=4)
#' dim(D)
#' coverage(D,4)
#' eCAN(4,8,2)   ## size larger than eCAN
#'
#' @importFrom sfsmisc digitsBase

#' @export
smcCA <- function(t, k, v, start0=TRUE, ...){
  Call <- sys.call()
  stopifnot(is.numeric(t), is.numeric(k), is.numeric(v))
  zeilen <- SMCcat[SMCcat$t>=t & SMCcat$k>=k & SMCcat$v==v,]
  ## lower strength is at the top,
  ## higher strength at the bottom
  zrows <- nrow(zeilen)
  if (zrows == 0) stop("The request cannot be served by function smcCA.")
  zeile <- zeilen[zrows + 1 - which.min(rev(zeilen$N)), ,drop=FALSE]
  vc <- as.character(zeile$v); tc <- as.character(zeile$t); kc <- as.character(zeile$k)
  aus <- SMC(SMC_CPHFs[[tc]][[vc]][[kc]], zeile$v, zeile$t)
  if (!start0) aus <- aus + 1
  aus[,1:k]
  class(aus) <- c("ca", class(aus))
  attr(aus, "Call") <- Call
  attr(aus, "origin") <- paste0("Sherwood, Martirosyan, Colbourn (2006), strength t=", t, ", ")
  attr(aus, "t") <- t
  aus
}

#' @export
SMC <- function(cphf, q, t, ...){
  # cphf a strength t SCPHF
  # in q^t levels
  ## matrix with the beta vectors in the order 0 top
  ## to t-1 bottom

  if (q %in% c(8,9)) gf <- mygf(q, 11) else gf <- lhs::create_galois_field(q)
  ## for 9, the primitive would be 17, but that does not make a difference
  ## the cphf for q=8 does not work with the field from lhs
  ##  but needs a different isomorphic variant of the gf

  ## coefficients for low powers at the top

  ## lspace provides the h tuples for the CA columns
  lspace <- sfsmisc::digitsBase(0:(q^(t-1)-1), q, ndigits=t-1)[(t-1):1,]
  ## betaspace contains the beta coefficients

  betaspace <- sfsmisc::digitsBase(0:(q^t-1), q, ndigits=t)[t:1,]

  CA <- gf_matmult(t(betaspace), rbind(1, lspace), gf)
  DHF2CA(cphf, CA)
}

#' @export
N_smcCA <- function(t, k, v){
  stopifnot(is.numeric(t), is.numeric(k), is.numeric(v))
  zeilen <- SMCcat[SMCcat$t>=t & SMCcat$k>=k & SMCcat$v==v,]
  ## lower strength is at the top,
  ## higher strength at the bottom
  zrows <- nrow(zeilen)
  if (zrows == 0) return(NA)
  zeile <- zeilen[zrows + 1 - which.min(rev(zeilen$N)), ,drop=FALSE]
  min(zeilen$N)
}

#' @export
k_smcCA  <- function(t, N, v){
  hilf <- SMCcat[SMCcat$t==t & SMCcat$N<=N &
                               SMCcat$v==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  max(hilf[,"k"])
}

SMC9 <- function(cphf, q, t, primitive, ...){
  ## experiment with the unsolved cphf
  ## no longer needed, keep for similar cases

  # cphf a strength t SCPHF
  # in q^t levels
  ## matrix with the beta vectors in the order 0 top
  ## to t-1 bottom

  if (!q==9) stop("for q=9 only")
  gf <- mygf(q, primitive)
  print(switch(as.character(primitive), '14'="default GF",
               '10'="x^2 + 1", '17'="x^2 + 2x + 2", "???"))

  ## lspace provides the h tuples for the CA columns
  lspace <- sfsmisc::digitsBase(0:(q^(t-1)-1), q, ndigits=t-1)[(t-1):1,]
  ## betaspace contains the beta coefficients

  betaspace <- sfsmisc::digitsBase(0:(q^t-1), q, ndigits=t)[t:1,]

  CA <- gf_matmult(t(betaspace), rbind(1, lspace), gf)
  DHF2CA(cphf, CA)
}
