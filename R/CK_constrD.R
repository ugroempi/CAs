#' Function to augment a design with Chateauneuf-Kreher construction D
#'
#' A strength 3 CA in k (v-1)-level columns and N runs is extended
#' to k v-level columns in N + k*M + k*(v-1) runs using a strength 2 CA in
#' M runs for k-1 (v-1)-level columns (construction D)
#'
#' @rdname CK_constrD
#'
#' @aliases CK_constrDCA
#' @aliases CK_constrD
#' @aliases N_CK_constrDCA
#'
#' @usage CK_constrDCA(k, v, ...)
#' @usage CK_constrD(D3, D2=NULL, check=FALSE, start0=TRUE, ...)
#' @usage N_CK_constrDCA(t=3, k, v, ...)
#'
#' @param k integer: the target number of columns
#' @param v integer: the number of levels
#' @param t integer: the target strength (always 3)
#' @param D3 a strength 3 CA with k columns in v-1 levels, coded as integers
#'           0, ..., v-2 or 1,...,v-1
#' @param D2 a strength 2 CA with k-1 columns in v-1 levels, coded as integers
#'           0, ..., v-2 or 1,...,v-1, same coding as \code{D3};\cr
#'           \code{NULL} implies automatic creation with function \code{\link{bestCA}},
#'           which uses the known optimal function \code{\link{KSK}} for \code{v}=2,
#'           and the best implemented construction for larger \code{v}.
#' @param check a logical: if TRUE, the strength requirements are verified
#' (may take very long for large arrays, therefore defaults to FALSE)
#' @param start0 logical: Do integer values start with 0 ?  (\code{D3} and \code{D2}
#'           must be compatible)
#' @param ... further arguments to function \code{\link{coverage}} (\code{CK_constrD})
#'       or \code{\link{bestN}} (function \code{N_CK_constrDCA})
#'
#' @details
#' \code{CK_constrDCA} is the interface function that uses the workhorse
#' function \code{CK_constrD}.
#'
#' @returns Function \code{CK_constrDCA} returns a matrix of class \code{ca}
#' with attributes, which is a strength 3 covering array.\cr
#' Function \code{CK_constrD} returns a matrix without class and attributes,\cr
#' function \code{N_CK_constrDCA} returns the run size of the CA.
#'
#' Note that the run size may depend on the availability of an internet connection
#' for cases for which the best strength 3 CA for \code{k} columns or the
#' best strength 2 CA for \code{k-1} columns is
#' in the Dwyer data base or the NIST library only
#' (these might exist - has not been checked).
#'
#' @references Chateauneuf and Kreher (2002)
#'
#' @examples
#' N_CK_constrDCA(3,6,6)
#' bestN(3,6,6) ## not competitive
#' dim(CK_constrD(SCA_Busht(5,3)))
#' dim(CK_constrDCA(6,6)) ## always strength 3
#'
#' @export
CK_constrD <- function(D3, D2=NULL, check=FALSE, start0=TRUE, ...){
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
    stopifnot(ncol(D2)>=k-1)
    ll <- levels.no(D2)
    stopifnot(all(ll==v))
    if (check) stopifnot(all(coverage(D2,2)==1))
    if (start0) stopifnot(min(D2)==0) else stopifnot(min(D2)==1)
  }
  ## provide the best implemented CA for D2
  if (is.null(D2)) D2 <- bestCA(2, k-1, v) + as.numeric(!start0)
  if (start0){
    D3 <- D3 + 1
    D2 <- D2 + 1
  }
  M <- nrow(D2)
  V <- 1:v
  II <- matrix(c(rep(0,M), rep(D2[,1], k-1)), k*M, 1)
  for (i in 1:(k-2))
    II <- cbind(II, c(rep(D2[,i],i), rep(0,M), rep(D2[,i+1], k-i-1)))
  ## last column
  II <- cbind(II, c(rep(D2[,k-1],k-1),rep(0,M)))
  III <- matrix(c(V, rep(0, (k-1)*v)), k*v, 1)
  for (i in 2:k)
    III <- cbind(III, c(rep(0,(i-1)*v),V,rep(0,(k-i)*v)))
  aus <- rbind(D3, II, III)
  if (!start0) aus <- aus + 1
  return(aus)
}

#' @export
CK_constrDCA <- function(k, v, ...){
  ## t=3 always holds
  if (v==2) stop("augmentation with constrD does not work for v=2")
  Call <- sys.call()
  k3 <- k
  D3 <- bestCA(3, k3, v-1)
  D2 <- bestCA(2, k3-1, v-1)  ## strength 2
  aus <- CK_constrD(D3, D2)
  class(aus) <- c("ca", class(aus))
  attr(aus, "origin") <- "Construction D, Chateauneuf and Kreher 2002"
  attr(aus, "Call") <- Call
  attr(aus, "t") <- 3
  aus
}

#' @export
N_CK_constrDCA <- function(t=3, k, v, ...){
  if (!t==3) return(NA)
  if (v==2) return(NA)
  ## N + k*M + k*(v-1) runs
  unname(bestN(3, k, v-1) +
           k*bestN(2, k-1, v-1) +
           k*(v - 1))
}
