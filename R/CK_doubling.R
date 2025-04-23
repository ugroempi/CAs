#' #' Function to create a strength 3 CA with Chateauneuf-Kreher doubling
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
#' @param D3 a strength 3 CA, integers 0, ..., v-1 or 1,...,v
#' @param D2 a strength 2 CA, integers 0, ..., v-1 or 1,...,v, same coding as \code{D3}; for v=2, \code{NULL} implies automatic creation with function \code{\link{KSK}}.
#' @param check a logical: if TRUE, the strength requirements are verified
#' (may take very long for large arrays, therefore defaults to FALSE)
#' @param start0 logical: Do integer values start with 0 ?
#' @param ... further arguments to function \code{\link{coverage}}
#'
#' @details
#' Chateauneuf and Kreher's (2002) doubling contributes a few best-known
#' strength 3 CAs to the Colbourn tables; in particular, where v is a prime
#' or prime power and thus fulfills the conditions for the Bose
#' construction of a strength 2 OA and the Bush
#' construction of a strength 3 OA (as implemented in R package lhs),
#' the designs are at present (13 Feb 2025) best-known for v>=5
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
#' ## Thm 4.5
#' D3 <- cyc(37,4,type="4a")
#' dim(D3)
#' ## do not have good strength 2 4-level CA of strength 2
#'
#' D3 <- DoE.base::L64.4.6 - 1 ## 6 columns are not created from lhs::createBush
#' D2 <- lhs::createBose(5,6,bRandom=FALSE)
#' ## fix 1 symbols (Meagher-Stevens) is best
#' ##    according to the Colbourn tables, with 19 runs
#' ## applying fuse to the above D2 yields 25-2=23 runs
#' D2 <- fuse(D2, 5, 4)  ## 64 + 3*23, 12 more than necessary
#' dim(D2)
#' E <- CK_doubling(D3, D2)
#' dim(E)
#' coverage(E, 3)
#' eCAN(3, 12, 4)
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
  if (is.null(D2)) if (!v==2) stop("D2=NULL requires v=2")
  if (v==2 && is.null(D2))
    D2 <- KSK(k=k) + as.numeric(!start0)
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

