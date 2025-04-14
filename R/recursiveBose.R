#' function to increase numbers of columns for strength 2 CA in q levels (q prime power) from k to kq+1
#'
#' applies construction principle of Hartman (2005) Thm 7.1 and Corollary 7.2
#'
#' @rdname recursiveBose
#'
#' @aliases recursiveBose
#'
#' @usage recursiveBose(q, A=NULL, d=NULL, ...)
#'
#' @param q a prime or prime power number of levels
#' @param A \code{NULL}, or a uniform CA of strength 2 in N rows and k columns with integer-valued levels (0 to q-1 or 1 to q). \code{A} takes precedence over \code{d}.
#' @param d \code{NULL}, or a small positive integer; 1 returns the Bose OA (q^2 x (q+1)), d>1 the processed Bose OA in d*(q^2) - (d-1)*q runs and q^d+q^(d-1)+...+q+1 columns.\cr If \code{A} is specified, \code{d} is ignored.
#' @param ... currently not used
#'
#' @returns a uniform strength 2 CA with \code{q} level columns.\cr
#' If neither \code{A} nor \code{d} are specified,
#' the result has \code{2q^2 - q} rows and \code{q^2 + q + 1} columns,
#' levels are coded from 0 to \code{q - 1}.\cr
#' If \code{A} is specified, the result has \code{N + q^2 - q} rows and \code{q*k + q + 1} columns,
#' in the same coding as \code{A}.\cr
#' If \code{d} is specified and \code{A=NULL},
#' the result has \code{d*q^2 - (d-1)*q} rows and \code{sum_i=0^d q^i} columns,
#' coded from 0 to \code{q}-1.\cr
#'
#' @section Details:
#' The function uses a method presented by Hartman (2005) (Theorem 7.1 and Corollary 7.2).
#' For \code{A=NULL}, the function uses `lhs::createBose(q, q+1, bRandom=TRUE)` for creating A
#' The user can specify any uniform strength 2 CA for \code{A}; the extension to more columns and runs will always use the Bose OA without the first q rows, and thus always add \code{q^2-q} rows.
#'
#' @examples
#' dim(D <- recursiveBose(3))
#' eCAN(2,13,3)  ## optimal
#' dim(recursiveBose(3, d=2))
#' dim(D2 <- recursiveBose(3, D))
#' eCAN(2,40,3)  ## reasonably close to optimal
#' dim(recursiveBose(3, d=3))
#' dim(D3 <- recursiveBose(3, D2))
#' eCAN(2,121,3)  ## still reasonably close to optimal
#' dim(recursiveBose(3,d=4))
#'
#' ## constructions with A=NULL are optimal
#' for (q in c(4,5,7,9,11,13,16)){
#'   D <- recursiveBose(q)
#'   cat(paste0("q=", q, "\n"))
#'   cat(paste0("constructed: ",
#'      paste(paste0(c("N=","ncolumns="),dim(D)), collapse=", "), "\n"))
#'   Nopt <- eCAN(2,ncol(D),q)[1]   ## optimal
#'   cat(paste0("optimal for ", ncol(D), " columns: ", Nopt, " runs\n"))
#' }
#'
#' #' ## constructions of first recursion are optimal or reasonably close
#' for (q in c(4,5,7,9,11,13,16)){
#'   D <- recursiveBose(q)
#'   D <- recursiveBose(q, D)
#'   cat(paste0("q=", q, "\n"))
#'   cat(paste0("constructed: ",
#'      paste(paste0(c("N=","ncolumns="),dim(D)), collapse=", "), "\n"))
#'   Nopt <- eCAN(2,ncol(D),q)[1]   ## optimal
#'   cat(paste0("optimal for ", ncol(D), " columns: ", Nopt, " runs\n"))
#' }
#'
#' ## a general A, an OA in 11 5-level columns at strength 2
#' ## coded starting with 1
#' A <- matrix(as.numeric(as.matrix(
#'      DoE.base::oa.design(nlevels=rep(5,11), randomize=FALSE))),
#'    nrow=50)
#' head(A)
#' ## number of columns increased from 11 to 56
#' ## number of rows increased by 5^2 - 5 = 20
#' ## still strength 2
#' dim(D <- recursiveBose(5, A))
#' coverage(D, 2, start0=FALSE)
#' ## 25% more runs than the best CA
#' eCAN(2,56,5)
#'
#' ## recursive with d=1
#' recursiveBose(3,d=1)
#'
#' @export
recursiveBose <- function(q, A=NULL, d=NULL, ...){
  start0 <- TRUE ## default
  if (!is.null(A)) d <- NULL
  if (!is.null(d)){
    stopifnot(is.numeric(d))
    stopifnot(d>0)
    stopifnot(d %% 1 == 0)
  }
  ## t=2
  ## q^2+q+1 columns
  ## 2q^2-q runs
  ## Hartman 2005, Thm 7.1 and Corollary 7.2
  C <- lhs::createBose(q, q+1, bRandom=FALSE)  ## checks for prime
  if (is.null(A)) A <- C else{
    stopifnot(is.matrix(A))
    stopifnot(is.numeric(A))
    if (min(A)==1 && max(A)==q) {
      start0 <- FALSE
      A <- A - 1
    }
    stopifnot(all(A %in% 0:(q-1)))
    if(length(unique(levels.no(A)))>1) stop("All columns of A must have the same number of levels.")
    stopifnot(all(coverage(A,2)==1))
  }
  if (is.null(d)){
    k <- ncol(A)
    R <- C[-(1:q), ]  ## omit the rows that have constants in columns 2:(q+1)
    # Williams: N <- rep(1:(q-1), each=q); this is column 1 of R
    top <- cbind(0, A[,rep(1:k, each=q)])
    bottom <- R[,c(1,rep(2:(q+1), k))]
    D <- rbind(top,bottom) + ifelse(start0, 0, 1)
  }else{
    ## recursive construction
    if (d==1) return(C)
    stopifnot(identical(A, C))
    D <- C
    for (i in 2:d){
      D <- recursiveBose(q, A)
      A <- D
    }
  }
  D
}
