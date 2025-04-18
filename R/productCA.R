#' Product constructions
#'
#' Functions to construct a CA from a product construction of two ingoing CAs
#'
#' @rdname productCA
#'
#' @aliases productCA1
#' @aliases productCA2
#'
#' @usage productCA1(D1, D2, check=TRUE, ...)
#' @usage productCA2(A, B, t, check=FALSE, ...)
#'
#' @param D1 an N x k CA of strength 2 with v levels
#' @param D2 an N x l CA of strength 2 with v levels
#' @param A  an N x k CA of strength \code{t} with v levels,
#' @param B  an M x k CA of strength \code{t} with w levels,
#' @param t strength of \code{A} and \code{B}
#' @param check logical, if TRUE, checks required strength of ingoing CAs (may substantially increase run time for larger CAs, especially for \code{t}>2)
#' @param ... currently not used
#'
#' @section Details:
#' \bold{I don't like the function names, they will likely change.}
#'
#' Function \code{productCA1} yields a CA(N+M, 2, kl, v) from a product construction, as, e.g., described in Torres-Jimènez et al. (2019).\cr
#' Function \code{productCA2} yields a CA(N*M, t, k, v*w) by combining the corresponding columns of the two CAs into \code{v*w} levels, according to Theorem 3.2 of Chateauneuf, Colbourn and Kreher (.\cr
#'
#' @returns Function \code{productCA1} returns a CA(N+M, 2, k*l, v) (N+M x k*l matrix).\cr
#' Function \code{productCA2} yields a CA(N*M, t, k, v*w) (N*M x k matrix).\cr
#' The returned matrices have their levels coded as the ingoing matrices, i.e., integers starting with 0 or 1.
#'
#' @references Torres-Jimenez et al. (2019), Chateauneuf and Kreher (2002)
#'
#' @examples
#' # product CA (N+M x k*l matrix)
#' A <- cyc(19,2)
#' B <- KSK(k=12)
#' dim(A); dim(B)
#' E <- productCA1(A, B)
#' dim(E)
#' coverage(E, 2)
#' ## more than twice as large as optimal
#' eCAN(2,228,2)
#'
#' ## two CAs with four columns each and strength 3
#' A <- lhs::createBush(3, 4, bRandom=FALSE)  ## strength 3, 4 columns in 3 levels
#' dim(A)
#' B <- lhs::createBush(4, 4, bRandom=FALSE)   ## strength 3, 4 columns in 4 levels
#' dim(B)
#' ## product array with 3*4=12 levels and 27*64=1728 runs
#' E <- productCA2(A, B, 3)
#' dim(E)
#' tail(E)
#' coverage(E,3)
#' eCAN(3,4,12)   ## the CA is optimal (as of April 2025)
#'
#'
#' @export
productCA1 <- function(D1, D2, check=TRUE, ...){
   stopifnot(is.matrix(D1)); stopifnot(is.matrix(D2))
   levs1 <- unique(levels.no(D1))
   levs2 <- unique(levels.no(D2))
   stopifnot(length(levs1)==1);stopifnot(length(levs2)==1)
   stopifnot(levs1==levs2)
   stopifnot(all(unique(c(D1))==unique(c(D2))))
   if(check){
     stopifnot(all(coverage(D1, 2)==1))
     stopifnot(all(coverage(D2, 2)==1))
   }
   k <- ncol(D1); l <- ncol(D2)
   aus <- rbind(D1[,rep(1:k, times=l)],
                D2[,rep(1:l, each=k)])
   return(aus[!duplicated(aus),])
}

#' @export
productCA2 <- function(A, B, t, check=FALSE, ...){
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
   aus
}
