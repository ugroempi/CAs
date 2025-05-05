#' Cross two strength t CAs with k columns to yield strength t in v*w levels
#'
#' Function to construct a strength t CA in k columns at v*w levels from a strength t CA in k columns at v levels and a strength t CA in k columns at w levels
#'
#' @rdname crossCAs
#'
#' @aliases crossCAs
#'
#' @usage crossCAs(A, B, t, check=FALSE, ...)
#'
#' @param A  an N x k CA of strength \code{t} with v levels,
#' @param B  an M x k CA of strength \code{t} with w levels,
#' @param t strength of \code{A} and \code{B}
#' @param check logical; if TRUE, checks required strength of ingoing CAs (may substantially increase run time for larger CAs, especially for \code{t}>2)
#' @param ... currently not used
#'
#' @section Details:
#' Function \code{crossCAs} yields a CA(N*M, t, k, v*w) by combining the corresponding columns of the two CAs
#' into \code{v*w} levels, according to, e.g., Theorem 2.8 of Zhang et al. (2014).
#'
#' @returns Function \code{crossCAs} yields a CA(N*M, t, k, v*w) (N*M x k matrix of class \code{ca}).\cr
#' The returned matrix has its levels coded like the ingoing matrices, i.e., integer starting with 0 or 1.
#'
#' @references Zhang et al. (2014)
#'
#' @examples
#' ################################################
#' ## crossCAs
#' ################################################
#' ## two CAs with four columns each and strength 3
#' A <- lhs::createBush(3, 4, bRandom=FALSE)  ## strength 3, 4 columns in 3 levels
#' dim(A)
#' B <- lhs::createBush(4, 4, bRandom=FALSE)   ## strength 3, 4 columns in 4 levels
#' dim(B)
#' ## product array with 3*4=12 levels and 27*64=1728 runs
#' E <- crossCAs(A, B, 3)
#' dim(E)
#' tail(E)
#' coverage(E, 3)
#' eCAN(3, 4, 12)      ## the CA is optimal (as of April 2025)
#' eCAK(3, 1728, 12)   ## the optimal strength 3 12-level CA accommodates two more columns,
#'                     ## and even in an OA (12 is not a prime power, thus 13 columns not possible)
#'

#' @export
crossCAs <- function(A, B, t, check=FALSE, ...){
  Call <- sys.call()
   stopifnot(is.matrix(A)); stopifnot(is.matrix(B))
   stopifnot(is.numeric(A)); stopifnot(is.numeric(B))
   stopifnot(min(A)>=0); stopifnot(min(B)>=0)
   stopifnot(min(A)%%1==0); stopifnot(min(B)%%1==0)
   start0 <- TRUE
   if (min(A)==1){
     if(!min(B)==1) stop("Levels of A and B must both start with 0 or both start with 1")
     start0 <- FALSE
   } else
   if(!min(B)==0) stop("Levels of A and B must both start with 0 or both start with 1")
   v <- unique(levels.no(A))
   w <- unique(levels.no(B))
   stopifnot(length(v)==1); stopifnot(length(w)==1)
   k <- ncol(A); l <- ncol(B); stopifnot(k==l)
   N <- nrow(A); M <- nrow(B)
   if(check){
     stopifnot(all(coverage(A, t)==1))
     stopifnot(all(coverage(B, t)==1))
   }
   aus <- matrix(NA, N*M, k)
   if (!start0) {
     A <- A - 1
     B <- B - 1
   }
   for (i in 1:k)
      aus[,i] <- rep(A[,i], each=M)*w + rep(B[,i], v)
   if (!start0) aus <- aus+1
   class(aus) <- c("ca", class(aus))
   attr(aus, "Call") <- Call
   aus
}
