#' Construction in Hartman Theorem 7.4
#'
#' The number of columns of an N x k design A can be squared; the design size is usually not optimal.
#'
#' @rdname Hartman74
#'
#' @aliases Hartman74
#' @aliases N_Hartman74
#' @aliases k_Hartman74
#' @aliases Turan
#'
#' @usage Hartman74(A, v, t, B=NULL, ...)
#' @usage N_Hartman74(N_A, k_A, v, t, ...)
#' @usage k_Hartman74(N_A, k_A, v, t, ...)
#' @usage Turan(t, v)
#'
#' @param A an N x k CA of strength \code{t} with columns in \code{v} levels each, denoted as 0 to v-1 or 1 to v
#' @param v the number of levels (or the number of parts for \code{Turan})
#' @param t the strength (or the number of vertices for \code{Turan})
#' @param B NULL, or a strength 2 CA in at least \code{Turan(t,v)} columns with levels from 1 to k (levels from 0 to (k-1) are acceptable too). If not specified by the user, the Bose construction from \code{\link[lhs]{createBose}} is used, which yields at most k+1 columns, for a prime power k. If k is not a prime power or k+1 columns are not enough, the user must provide a suitable \code{B} or an \code{A} with more columns (may not be easy, at least not for achieving acceptable run sizes).
#' @param N_A run size of the matrix \code{A} (called N above)
#' @param k_A number of columns of the matrix \code{A} (called k above)
#' @param ... currently not used
#'
#' @section Details:
#' The construction vertically stacks \code{Turan(t,v)} copies of column selections from \code{A},\cr
#' where the ith copy is obtained by picking the \code{A} columns according to the elements of the ith column of \code{B}.\cr
#' Eventually, duplicate rows are removed.
#'
#' The construction requires that the number of columns of \code{B} is at least \code{Turan(t,v)}, which may be impossible if \code{A} has too few columns. Sometimes, it can be possible to increase the number of columns for A (see examples).
#'
#' @returns function \code{Hartman74} returns a strength \code{t} CA with \code{k^2} columns (matrix of class \code{ca}),\cr
#' the identical functions \code{N_Hartman74} and \code{k_Hartman74} return a named vector with
#' the run size \code{N}, the number of columns \code{k}, the strength \code{t} and the number of levels \code{v}\cr
#' from using an \code{N_A} x \code{k_A} matrix \code{A} of strength \code{t} in \code{v} levels,\cr
#' assuming this matrix and a suitable B exists (warning, if this is not the case for the Bose
#' construction for prime power k),\cr
#' and function \code{Turan} returns the maximum number of edges in a \code{v}-partite graph
#' with \code{t} vertices, which is used in the construction.
#'
#' @references Hartman (2005)
#'
#' @examples
#' # Turan numbers (minimum number of columns needed for B)
#' Turanvec <- Vectorize(Turan)
#' TuranTable <- outer(1:13, 2:7, function(X, Y) Turanvec(X, Y))
#' dimnames(TuranTable) <- list(t=1:13, v=2:7)
#' TuranTable
#' ## the first seven rows remain constant for larger values of v
#'
#' A <- lhs::createBusht(4,5,3,bRandom=FALSE) ## 64 x 5
#' D <- Hartman74(A, v=4, t=3)
#' coverage(D, 3)
#' dim(D)
#' eCAN(3,25,4) ## D is far from optimal (roughly factor 1.6)
#'
#' A <- lhs::createBusht(3,4,4,bRandom=FALSE) ## 81 x 4
#' A5 <- derive(DoE.base::L243.3.6-1)         ## 81 x 5
#'      ## optimal CA for strength 4 in 5 3-level columns according to Colbourn tables
#' # D <- Hartman74(A, v=3, t=4) ## would yield an error, too few columns in Bose design
#' D <- Hartman74(A5, v=3, t=4)
#' # coverage(D, 4)  ## commented out because of check time
#' dim(D)
#' eCAN(4,25,3) ## D is far from optimal (roughly factor 1.44)
#'
#' @export
N_Hartman74 <- function(N_A, k_A, v, t, ...){
  ## N_A the run size of the ingoing array A
  ## v the number of levels
  ## t the strength of the ingoing array A
  ## existence of A is assumed
  stopifnot(is.numeric(N_A))
  stopifnot(N_A>0)
  stopifnot(N_A%%1==0)
  stopifnot(is.numeric(k_A))
  stopifnot(k_A>0)
  stopifnot(k_A%%1==0)
  stopifnot(is.numeric(v))
  stopifnot(v>0)
  stopifnot(v%%1==0)
  stopifnot(is.numeric(t))
  stopifnot(t>0)
  stopifnot(t%%1==0)
  if (k_A+1 < Turan(t,v)) warning("The returned size is likely not achievable.")
  c(N=Turan(t,v)*N_A, k=k_A^2, v=v, t=t)
}
#' @export
k_Hartman74 <- N_Hartman74

#' @export
Hartman74 <- function(A, v, t, B=NULL, ...){
  Call <- sys.call()
  ## A a matrix of strength t in k v level columns
  ## B a matrix of strength 2 in Turan(t,v)+1 k-level columns, or NULL
  stopifnot(v>0)
  stopifnot(v%%1==0)
  stopifnot(t>0)
  stopifnot(t%%1==0)
  stopifnot(is.matrix(A))
  stopifnot(is.numeric(A))
  if (!is.null(B)){
    stopifnot(is.matrix(B))
    stopifnot(is.numeric(B))
  }
  start0 <- TRUE
  if (min(A)==1){
    start0 <- FALSE
    A <- A - 1
  }
  stopifnot(all(A %in% 0:(v-1)))
  N_A <- nrow(A); k <- ncol(A)
  colBneeded <- Turan(t,v)+1

  if (is.null(B)){
    B <- try(lhs::createBose(k,k+1), silent = TRUE)
    if("try-error" %in% class(B)) stop("If k is not a prime power, B must be provided.")
    if (ncol(B) < colBneeded) {
      message("the default B has too few columns, and an adequately-sized B for this A is likely impossible")
      message("you can try to use an A with more columns, which will imply more columns in the result")
      stop("too few columns in B,\n at least ", colBneeded, " columns needed")
    }
  } else {
    if (min(B)==1) B <- B - 1
    stopifnot(all(B %in% 0:(v-1)))
    if (ncol(B)<colBneeded) stop("too few columns in B,\nat least ", colBneeded, " columns are needed")
  }
  B <- B + 1  ## because to be used as position index in A

  Cs <- lapply(1:colBneeded, function(obj){
    A[,B[,obj]]
  })
  aus <- do.call(rbind,c(as.list(0:(v-1)),Cs))
  aus <- aus[!duplicated(aus),]
  class(aus) <- c("ca", class(aus))
  attr(aus, "Call") <- Call
  if (start0) return(aus) else return(aus+1)
}

#' @export
Turan <- function(t, v){
  # calculates the maximum number of edges in a graph with t vertices and n parts
  # t number of vertices
  # v number of parts (called n in Hartman paper)
  if (v >= t) return(t*(t-1)/2)  ## complete graph
  ## v < t, i.e. t/v > 1
  ## na parts of size a and nb parts of size b
  a <- floor(t/v); b <- a + 1; nb <- t - a*v; na <- v - nb
  ## subtract number of edges within parts
  t*(t-1)/2 - na*a*(a-1)/2 - nb*b*(b-1)/2
}
