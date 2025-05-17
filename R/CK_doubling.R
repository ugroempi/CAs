#' Function to create a strength 3 CA with Chateauneuf-Kreher doubling
#'
#' A strength 3 CA in k v-level columns and N runs is extented
#' to 2k v-level columns in N + M*(v-1) runs
#' by adding cyclic permutations of a strength 2 CA in v levels
#' and M runs (Thms 4.1 and 4.5)
#'
#' @rdname CK_doubling
#'
#' @aliases CK_doubling
#'
#' @usage CK_doubling(D3, D2=NULL, check=FALSE, start0=TRUE, ...)
#'
#' @param D3 a strength 3 CA with k columns in v levels, coded as integers
#'           0, ..., v-1 or 1,...,v
#' @param D2 a strength 2 CA with k columns in v levels, coded as integers
#'           0, ..., v-1 or 1,...,v, same coding as \code{D3};\cr
#'           for v=2 or v=3, \code{NULL} implies automatic creation with function
#'           \code{\link{KSK}} (v=2) or function \code{\link{CAEX}} (v=3).
#' @param check a logical: if TRUE, the strength requirements are verified
#' (may take very long for large arrays, therefore defaults to FALSE)
#' @param start0 logical: Do integer values start with 0 ?  (\code{D3} and \code{D2}
#'           must be compatible)
#' @param ... further arguments to function \code{\link{coverage}}
#'
#' @details
#' Chateauneuf and Kreher's (2002) doubling contributes a few best-known
#' strength 3 CAs to the Colbourn tables; in particular, where v is a prime
#' or prime power and thus fulfills the conditions for the Bose
#' construction of a strength 2 OA and the Bush
#' construction of a strength 3 OA (as implemented in R package lhs),
#' the designs are at present (13 Feb 2025) best-known for v >= 5
#'
#' @references Chateauneuf and Kreher (2002)
#'
#' @examples
#' L8.2.4 <- cbind(rep(1:2,each=4), rep(1:2, each=2,times=2), rep(1:2,4), c(1,2,2,1,2,1,1,2))
#' coverage(L8.2.4,3,start0=FALSE)
#' E <- CK_doubling(L8.2.4, start0=FALSE)
#' dim(E)
#' coverage(E-1, 3)  ## successful
#' eCAN(3,8,2)   ## not optimal, but close
#'
#' D3 <- paleyHad(11)
#' coverage(D3, 3)  # suitable
#' E <- CK_doubling(D3)
#' dim(E)
#' coverage(E, 3)  ## successful
#' eCAN(3,22,2)   ## optimal
#'
#' ## Thm 4.5, 4 levels
#' D3 <- DoE.base::L64.4.6 - 1 ## 6 columns are not created from lhs::createBush
#' D2 <- CS_MS(6, 4)           ## 19 runs, optimal
#' eCAN(2, 6, 4)
#' E <- CK_doubling(D3, D2)
#' coverage(E, 3)   ## has worked
#' dim(E)           ## size almost optimal
#' eCAN(3, 12, 4)
#'
#' ## 3 levels, with default CAEX-design
#' ## D3 is the current best-known CA(42, 3, 8, 3) from a cross-summing construction
#' ##     by Colbourn et al. (2010, CKRS)
#' eCAN(3, 8, 3)
#'
#' CA42.3.8.3 <- rbind(
#' c(0,0,0,0,0,0,0,0), c(1,1,1,1,1,1,1,1), c(2,2,2,2,2,2,2,2), c(1,2,2,0,1,1,0,1),
#' c(1,0,0,2,1,2,2,2), c(1,0,1,1,0,2,0,1), c(2,2,2,2,0,0,0,0), c(0,2,1,1,1,0,0,2),
#' c(0,2,1,0,2,2,1,1), c(0,2,0,2,0,2,1,2), c(0,1,0,1,2,0,1,0), c(2,0,1,2,0,0,1,1),
#' c(1,1,0,2,2,0,0,1), c(2,2,1,2,1,2,1,0), c(1,1,2,0,0,2,1,0), c(2,0,1,1,0,1,2,0),
#' c(2,0,0,1,2,2,1,2), c(2,1,1,0,2,0,0,0), c(2,1,2,2,1,0,1,2), c(0,0,0,0,1,1,1,1),
#' c(1,0,2,2,2,1,1,0), c(1,0,2,1,0,0,2,2), c(1,0,1,0,1,0,2,0), c(1,2,1,2,0,1,2,1),
#' c(0,1,2,0,1,1,2,2), c(2,2,1,0,0,1,1,2), c(0,0,2,0,2,0,2,1), c(2,2,0,1,1,0,2,1),
#' c(0,1,2,2,1,2,0,1), c(0,1,1,2,0,0,2,0), c(0,2,2,1,0,1,1,1), c(0,2,0,0,2,1,2,0),
#' c(1,1,1,1,2,2,2,2), c(2,1,0,0,0,2,2,1), c(2,1,0,2,1,1,0,0), c(2,1,2,1,2,1,0,1),
#' c(2,0,2,0,1,2,0,2), c(1,2,0,1,2,2,0,0), c(0,0,2,1,1,2,2,0), c(1,1,0,1,0,1,0,2),
#' c(0,0,1,2,2,1,0,2), c(1,2,0,0,2,0,1,2)
#' )
#'
#' ## D2=NULL implies use of CAEX(k=8), which has 13 runs
#' E <- CK_doubling(CA42.3.8.3)
#' dim(E)    ## 42 + 2*13 = 68 runs with 16 columns
#' # coverage(E, 3)  ## commented out because of example run time
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
  if (is.null(D2)) if (!v %in% c(2,3)) stop("D2=NULL requires v=2 or v=3")
  if (is.null(D2))
    if (v==2)
      D2 <- KSK(k=k) + as.numeric(!start0)
    else  ## v==3
      D2 <- CAEX(k=k) + as.numeric(!start0)
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

