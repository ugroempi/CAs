#' Cross-sum two codes with k columns in v levels
#'
#' Function to construct cross-sum two codes in k columns at v levels, according to CKRS (2010)
#'
#' @rdname crossSum
#'
#' @aliases crossSum
#' @aliases directSum
#'
#' @usage crossSum(A, B, ...)
#' @usage directSum(A, B, ...)
#'
#' @param A  an N x k code with (up to) \code{v} non-negative levels.
#' @param B  an M x l code with (up to) \code{v} non-negative levels;\cr
#'           for function \code{crossSum}, \code{l=k} is needed.
#' @param ... currently not used
#'
#' @section Details:
#' The convenience function \code{directSum} yields
#' \code{cbind(kronecker(A, matrix(1,M,1)), kronecker(matrix(1,N,1), B)}, i.e.,
#' it concatenates all pairs of rows of matrices \code{A} and \code{B}.
#'
#' Function \code{crossSum} yields an \code{N*M x k} matrix by adding all pairs of
#' rows from \code{A} and \code{B} (each row of \code{A} to all columns of \code{B}).
#' According to CKRS (2010), it is most promising to choose \code{A} as a repetition code (RC) (e.g.,
#' with rows 0000 1111 2222, for with k=4 and v=3), an RC extended by a
#' few constant columns (ERC, e.g., 000000 111100 222200),
#' a direct sum of repetition codes (DRC, e.g., rows
#' 0000000000000000 1111111111111111 2222222222222222 1111111111111111
#' 2222222222222222 0000000000000000 2222222222222222 0000000000000000
#' 1111111111111111 for the above RC with itself),
#' or a DRC extended by a few constant columns EDRC.
#'
#' @returns \code{directSum} returns an \code{N*M x (k+l)} matrix.\cr
#' \code{crossSum} returns an \code{N*M x k} matrix, which has its levels
#' coded like the ingoing matrices, i.e., integer starting with 0 or 1.
#'
#' @references Colbourn, Keri, Rivas Soriani and Schlage-Puchta (2010, CKRS)
#'
#' @examples
#' ################################################
#' ## crossSum
#' ################################################
#' ## two codes with fourteen columns each
#' ## A a repetition code
#' A <- matrix(0:2, 3, 14)
#' dim(A)
#' ## B from Fig. 1 of CKRS
#' B <- rbind(
#' c(0,0,0,0,0,0,0,0,0,0,0,0,0,0), c(0,2,2,2,1,1,2,0,2,0,0,2,1,1),
#' c(2,0,2,0,1,2,0,1,0,1,1,1,0,2), c(1,0,1,2,2,2,1,1,0,0,0,2,0,0),
#' c(2,0,0,2,0,0,1,1,2,1,2,1,2,1), c(0,0,1,1,0,2,1,2,2,2,1,0,1,2),
#' c(1,2,0,0,0,2,1,2,1,0,0,1,2,2), c(0,0,1,2,1,0,0,2,1,0,1,2,2,1),
#' c(2,2,1,0,2,0,1,0,0,2,1,1,2,1), c(1,2,0,2,1,0,2,1,1,2,1,0,0,2),
#' c(0,2,2,0,1,0,1,1,2,0,0,2,1,2), c(0,0,1,2,2,0,2,0,1,1,0,1,1,2),
#' c(2,1,1,1,2,0,0,1,2,0,0,0,2,2), c(0,1,0,0,2,2,2,1,2,0,1,1,0,1),
#' c(0,2,1,0,1,2,2,1,1,1,2,0,2,0), c(0,2,0,2,2,2,0,0,2,1,1,2,2,2),
#' c(0,2,1,2,0,0,1,0,2,2,2,2,0,0)
#' )
#'
#' dim(B)
#' ## cross-sum array with 3*17=51 runs
#' E <- crossSum(A, B)
#' dim(E)
#' coverage(E, 3)
#' eCAN(3, 14, 3)      ## the CA is not optimal
#' eCAK(3, 51, 3)      ## the optimal strength 3 3-level CA
#'                     ## in 51 runs accommodates two more columns

#' @export
crossSum <- function(A, B, ...){
   stopifnot(is.matrix(A)); stopifnot(is.matrix(B))
   stopifnot(is.numeric(A)); stopifnot(is.numeric(B))
   stopifnot(min(A)>=0); stopifnot(min(B)>=0)
   stopifnot(min(A)%%1==0); stopifnot(min(B)%%1==0)

   ## set start0
   start0 <- TRUE
   if (min(A)==1){
     if(!min(B)==1) stop("Levels of A and B must both start with 0 or both start with 1")
     start0 <- FALSE
     A <- A-1; B <- B-1
   } else
   if(!min(B)==0) stop("Levels of A and B must both start with 0 or both start with 1")

   v <- unique(levels.no(B))
   stopifnot(length(v)==1); stopifnot(all(A %in% 0:(v-1)))
   k <- ncol(A)
   if (!(ncol(B)==k)) stop("A and B must have the same number of columns")

   N <- nrow(A); M <- nrow(B)
   aus <- matrix(NA, N*M, k)
   if (!start0) {
     A <- A - 1
     B <- B - 1
   }
   for (i in 1:N)
     for (j in 1:M)
      aus[(i-1)*M + j,] <- (A[i,] + B[j,])%%v
   if (!start0) aus <- aus+1
   class(aus) <- c("ca", class(aus))
   #attr(aus, "Call") <- Call
   aus
}

#' @export
directSum <- function(A, B, ...){
  stopifnot(is.matrix(A), is.matrix(B))
  stopifnot(is.numeric(A), is.numeric(B))
  N <- nrow(A); M <- nrow(B)
  k <- ncol(A); l <- ncol(B)
  cbind(kronecker(A, matrix(1,M,1)) , kronecker(matrix(1,N,1), B))
}
