#' Function to create a strength 3 CA with Chateauneuf-Kreher construction D
#'
#' A strength 3 CA in k v-level columns and N runs is extented
#' to 2k v-level columns in N + M*(v-1) runs
#' by adding cyclic permutations of a strength 2 CA in v levels
#' and M runs
#'
#' @rdname CK_constrD
#'
#' @aliases CK_constrD
#'
#' @usage CK_constrD(D3, D2=NULL, verify=FALSE, start0=TRUE, ...)
#'
#' @param D3 a strength 3 CA, integers 0, ..., v-1 or 1,...,v
#' @param D2 a strength 2 CA; if NULL, creates it by function \code{\link{KSK}}
#' (only for 2-level CAs); must have the same coding as \code{D3};
#' for prime or prime power v with up to v+1 columns, the user might use a Bose array
#' that can be created with function \code{\link[lhs]{createBose}}.
#' @param verify a logical: if TRUE, the strength requirements are verified
#' (may take very long for large arrays, therefore defaults to FALSE)
#' @param start0 logical: Do integer values start with 0 ?
#' @param ... further arguments to function \code{\link{coverage}}
#'
#' @details
#' Chateauneuf and Kreher's (2002) construction D does not necessarily yield smallest arrays.
#'
#' @examples
#' # create a CA(19,3,22,2) from a CA(12,3,11,2) and a CA(7,2,11,1)
#' A <- paleyHad(11)
#' B <- KSK(k=11)
#' CK_constrD(A, B)  ## explicitly specify B
#' CK_constrD(A)     ## this B is also used automatically (for v=2 only)
#' mistaken <- CK_constrD(A, B[-7,])  ## D2 does not have strength 2
#' coverage(mistaken, 3)      ## without verify, strength 3 is not guaranteed
#' # CK(A, B[-7,1], verify=TRUE)   ## would throw an informative error
#'
#' A <- lhs::createBush(3, 4, bRandom=FALSE)  ## strength 3 OA
#' B <- lhs::createBose(3, 4, bRandom=FALSE)   ## strength 2 OA
#' CK_constrD(A, B)                ## 45 runs, 8 3-level columns
#' CK_constrD(A+1, B+1, start0=FALSE)      ## levels starting with 1
#' eCAN(3, 8, 3)           ## 42 runs are smallest known size
#'
#' A <- lhs::createBush(5, 6, bRandom=FALSE)  ## strength 3 OA
#' B <- lhs::createBose(5, 6, bRandom=FALSE)   ## strength 2 OA
#' CK_constrD(A, B)                 ## 225 x 12
#' eCAN(3, 12, 5)           ## 225 runs are smallest known size
#' ## the analogous strategy also yields the best-known strength 3 CAs
#' ##    for v=7,8,9,11,13,16,17,19,23,25 with 2v+2 columns
#'
#' @export
CK_constrD <- function(D3, D2=NULL, verify=FALSE, start0=TRUE, ...){
  if (!is.matrix(D3)) D3 <- as.matrix(D3)
  if (!is.null(D2)){
    if (!is.matrix(D2)) D2 <- as.matrix(D2)
  }
  ll <- levels.no(D3)
  stopifnot(length(unique(ll))==1)
  v <- ll[1]
  if ((!min(D3)==0) && start0) stop("Did you forget to specify start0=FALSE?")
  if ((!start0) && !min(D3)==1) stop("With start0=FALSE, levels must be integers 1,...,v")
  if (!is.null(D2)) {
    ll <- levels.no(D2)
    stopifnot(all(ll==v))
  }
  stopifnot(length(setdiff(D2,D3))==0)
  if ((!v==2) && is.null(D2))
    stop("D2=NULL is permitted for v=2 levels only.")
  if (verify){
    stopifnot(coverage(D3,3, start0=start0)$min==1)
    if (!is.null(D2)) stopifnot(coverage(D2,2, start0=start0)$min==1)
  }
  if (is.null(D2)) D2 <- KSK(k=ncol(D3)) + (!start0)
  B <- do.call(rbind, rep(list(D2), v-1))
  Bs <- do.call(rbind, lapply(1:(v-1),
                              function(obj) (D2-(!start0)+obj)%%v + !start0))
  aus <- rbind(cbind(D3, D3),
        cbind(B, Bs))
  rownames(aus) <- NULL
  aus
}

