#' Function to create a strength 3 CA with Chateauneuf-Kreher doubling
#'
#' A strength 3 CA in j v-level columns and N runs is extended
#' to k=2j v-level columns in N + M*(v-1) runs
#' by adding cyclic permutations of a strength 2 CA in v levels
#' and M runs (Thms 4.1 and 4.5)
#'
#' @rdname CK_doubling
#'
#' @aliases CK_doublingCA
#' @aliases CK_doubling
#' @aliases N_CK_doublingCA
#'
#' @usage CK_doublingCA(k, v, ...)
#' @usage CK_doubling(D3, D2=NULL, check=FALSE, start0=TRUE, ...)
#' @usage N_CK_doublingCA(t=3, k, v, ...)
#'
#' @param k integer: the target number of columns
#' @param v integer: the number of levels
#' @param t integer: the target strength (always 3)
#' @param D3 a strength 3 CA with j columns in v levels, coded as integers
#'           0, ..., v-1 or 1,...,v
#' @param D2 a strength 2 CA with j columns in v levels, coded as integers
#'           0, ..., v-1 or 1,...,v, same coding as \code{D3};\cr
#'           \code{NULL} implies automatic creation with function \code{\link{bestCA}},
#'           which uses the known optimal function \code{\link{KSK}} for \code{v}=2,
#'           and the best implemented construction for larger \code{v}.
#' @param check a logical: if TRUE, the strength requirements are verified
#' (may take very long for large arrays, therefore defaults to FALSE)
#' @param start0 logical: Do integer values start with 0 ?  (\code{D3} and \code{D2}
#'           must be compatible)
#' @param ... further arguments to function \code{\link{coverage}} (\code{CK_doubling})
#'       or \code{\link{bestN}} (function \code{N_CK_doublingCA})
#'
#' @details
#' \code{CK_doublingCA} is the interface function that uses the workhorse
#' function \code{CK_doubling}.
#'
#' Chateauneuf and Kreher's (2002) doubling contributes a few best-known
#' strength 3 CAs to the Colbourn tables; in particular, where v is a prime
#' or prime power and thus fulfills the conditions for the Bose
#' construction of a strength 2 OA and the Bush
#' construction of a strength 3 OA (as implemented in \code{SCA_Busht}),
#' the designs are at present (13 Feb 2025) best-known for v >= 5.
#'
#' @returns Function \code{CK_doublingCA} returns a matrix of class \code{ca}
#' with attributes, which is a strength 3 covering array.\cr
#' Function \code{CK_doubling} returns a matrix without class and attributes,\cr
#' function \code{N_CK_doublingCA} returns the run size of the CA.
#'
#' Note that the run size may depend on the availability of an internet connection
#' for cases for which the best strength 3 CA for \code{ceiling(k/2)} columns is
#' in the Dwyer data base or the NIST library only
#' (these might exist - has not been checked).
#'
#' @references Chateauneuf and Kreher (2002)
#'
#' @examples
#' (E1 <- CK_doublingCA(8, 2))
#' coverage(E1, 3)  ## successful
#' eCAN(3,8,2)   ## not optimal, but close
#'
#' ## the above was achieved with these ingredients:
#' D3 <- bestCA(3, 4, 2)  ## 8x4, OA
#' dim(D3)
#' coverage(D3, 3)
#' D2 <- bestCA(2,4,2)   ## 5x4, KSK construction
#' dim(D2)
#' E2 <- CK_doubling(D3, D2)
#' dim(E2)
#' table(E1-E2) ## content is the same, attributes differ
#'
#' D3 <- paleyCA(3, 11)
#' coverage(D3, 3)  # suitable
#' E <- CK_doubling(D3)
#' dim(E)
#' coverage(E, 3)  ## successful
#' eCAN(3,22,2)   ## optimal
#'
#' E3 <- CK_doublingCA(12, 4)
#' dim(E3)
#' ## Thm 4.5, 4 levels
#' ## ingredients for E3:
#' D3 <- DoE.base::L64.4.6 - 1 ## 6 columns are not created from SCA_Busht
#'   ## bestCA(3, 6, 4) would yield a suitable CA from the DWYER database,
#'   ## provided that there is an internet connection
#' D2 <- bestCA(2, 6, 4)        ## 19 runs, optimal
#' eCAN(2, 6, 4)
#' ## E3 has 64 + (4-1)*19 = 121 rows
#' eCAN(3, 12, 4)    ## size almost optimal
#'
#' ## not all outcomes are so close to optimal
#' dim(CK_doublingCA(16, 3))
#' ## 42 + 2*13 = 68 runs
#' eCAN(3, 16, 3)  ## run size too large by factor 4/3 versus the current optimum 51
#'
#' @export
CK_doubling <- function(D3, D2=NULL, check=FALSE, start0=TRUE, ...){
  if (!is.matrix(D3)) D3 <- as.matrix(D3)
  ll <- levels.no(D3)
  stopifnot(length(unique(ll))==1)
  v <- ll[1]
  if (check) stopifnot(all(coverage(D3,3)==1))
  k <- ncol(D3)
  if (start0) stopifnot(min(D3)==0) else stopifnot(min(D3)==1)
  ## checks for non-NULL D2
  if (!is.null(D2)){
    if (!is.matrix(D2)) D2 <- as.matrix(D2)
    stopifnot(ncol(D2)>=k)
    ll <- levels.no(D2)
    stopifnot(all(ll==v))
    if (check) stopifnot(all(coverage(D2,2)==1))
    if (start0) stopifnot(min(D2)==0) else stopifnot(min(D2)==1)
  }
  ## provide the best implemented CA for D2
  if (is.null(D2)) D2 <- bestCA(2, k, v) + as.numeric(!start0)
  if (!start0){
    D3 <- D3 - 1
    D2 <- D2 - 1
  }
   list1 <- rep(list(D2), v-1)
   list2 <- lapply(1:(v-1), function(obj) (D2+obj) %% v)
    aus <- rbind(cbind(D3,D3),cbind(do.call(rbind, list1), do.call(rbind, list2)))
    if (!start0) aus <- aus + 1
    return(aus)
}

#' @export
CK_doublingCA <- function(k, v, ...){
  ## t=3 always holds
  Call <- sys.call()
  k3 <- ceiling(k/2)
  D3 <- bestCA(3, k3, v)
  D2 <- bestCA(2, k3, v)  ## strength 2
  aus <- CK_doubling(D3, D2)[,1:k]
  class(aus) <- c("ca", class(aus))
  attr(aus, "origin") <- "CK doubling, Chateauneuf and Kreher 2002"
  attr(aus, "Call") <- Call
  attr(aus, "t") <- 3
  aus
}

#' @export
N_CK_doublingCA <- function(t=3, k, v, ...){
  if (!t==3) return(NA)
  k3 <- ceiling(k/2)
  ## N + M*(v-1) runs
  unname(bestN(3, k3, v) +
    bestN(2, k3, v)*(v - 1))
}

# #' @export
# this is not up to use yet
# the problem is recursive use of CK_doubling
# maybe this must be directly addressed
k_CK_doublingCA <- function(t=3, N, v, ...){
  ## the best ingredient CAs are not expected to are themselves from
  ## CK_doubling
  if (!t==3) return(NA)
  hilfk <- ks(3, N, v); hilfk <- hilfk[-length(hilfk)] # omit eCAK
  kupper <- floor(max(hilfk)/2)
  ## N + M*(v-1) runs
  if (bestN(3, kupper, v) +
    bestN(2, kupper, v)*(v - 1) <= N) return(2*kupper)
  ## count down until feasible
  while(kupper<=1){
    kupper <- kupper - 1
    if (bestN(3, kupper, v) +
        bestN(2, kupper, v)*(v - 1) <= N) return(2*kupper)
  }
}
