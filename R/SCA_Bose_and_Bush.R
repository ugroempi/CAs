#' Functions for regular orthogonal arrays with prime power numbers of levels
#'
#' Functions for regular orthogonal arrays with prime power numbers of levels q.
#' SCA_Bose creates strength 2 OAs in q^2 runs with q+1 columns at q levels,
#' SCA_Busht creates strength t OAs (including t=2) in q^t runs with q+1 columns at q levels.
#'
#' @rdname SCA_Bose
#'
#' @aliases SCA_Bose
#' @aliases SCA_Busht
#'
#' @usage SCA_Bose(q)
#' @usage SCA_Busht(q, t)
#'
#' @param q integer: a prime power
#' @param t integer: the requested strength
#'
#' @section Details:
#' Function \code{SCA_Bose} obtains the strength 2 Bose constructed CA (q+1 columns in q levels, q^2 runs, q a prime power)
#' and rearranges it into SCA shape. The array actually is a 2PCA, in the sense that its bottom left \code{(q^2-q)xq} block is itself a PCA
#' with width of the left parts equal to \code{q} (and thus width of the right-hand side equal to zero). Every consecutive \code{qxq} block
#' of the matrix consists of columns whose elements are a permutation of the number 0 to \code{q-1}.
#' This structure of the matrix is exploited by function \code{\link{productPCA}}.
#'
#' Function \code{SCA_Busht} obtains the strength t Bush constructed CA (q+1 columns in q levels, q^t runs, q a prime power;
#' exception: for \code{t=3} and \code{q} a power of 2, $q+2$ columns are achieved.\cr
#' In all cases, the array is arranged into SCA shape with \code{k1=q}.
#' For \code{t=2}, it equals the array created by \code{SCA_Bose}.
#' This array is also a 2PCA; however, \code{productPCA} is a strength 2 construction only.
#'
#' Function \code{SCA_Bose} uses function \code{\link[lhs]{createBose}},
#' function \code{SCA_Busht} uses function \code{\link[lhs]{createBusht}},
#' except for \code{t=3} with \code{q} a power of 2, where it uses the internal
#' function \code{createBush3PowerOf2} of this package, which is based
#' on Sherwood et al. (2006).
#'
#' @returns
#' Both functions yield a matrix of class \code{ca} which is an orthogonal array in
#' \code{q+1} columns with \code{q} levels (exception: \code{q+2} columns for \code{t=3} and \code{q} a power of 2);
#' the number of rows is
#' \code{q^2} for \code{SCA_Bose} and \code{q^t} for \code{SCA_Busht}.
#' The arrays are rearranged such that their first \code{q} columns
#' have their first \code{q} rows constant,
#' which is beneficial for some CA constructions, e.g., for
#' \code{\link{powerCA}}.
#'
#' @references Bose (1938), Bush (1952), Sherwood, Martirosyan and Colbourn (2006)
#'
#' @examples
#' A <- SCA_Bose(4)  ## Bose array, arranged as SCA
#' B <- SCA_Busht(4,2)
#' all.equal(A, B)
#' ## arrays the same, difference only in attributes
#'
#' dim(SCA_Busht(4, 3))  ## q+2 columns for t=3 and q=2^s
#' dim(SCA_Busht(4, 4))  ## q+1 columns otherwise
#' dim(SCA_Busht(5, 3))  ## q+1 columns otherwise
#'

#' @export
SCA_Bose <- function(q){
  Call <- sys.call()
  aus <- lhs::createBose(q, q+1, bRandom=FALSE)[,c(2:(q+1),1)]
  class(aus) <- c("ca", class(aus))
  attr(aus, "t") <- 2
  attr(aus, "PCAstatus") <- list(type="SCA", k1=q, k2=1)
  attr(aus, "Call") <- Call
  attr(aus, "origin") <- "Bose construction"
  aus
}

#' @export
SCA_Busht <- function(q, t){
  Call <- sys.call()
  if (t==3 && ((q/2)%%1)==0) aus <- createBush3PowerOf2(q) else
  aus <- lhs::createBusht(q, q+1, t, bRandom=FALSE)[,c(2:(q+1),1)]
  k <- ncol(aus)
  class(aus) <- c("ca", class(aus))
  attr(aus, "Call") <- Call
  attr(aus, "t") <- t
  attr(aus, "PCAstatus") <- list(type="SCA", k1=q, k2=k-q)
  attr(aus, "origin") <- "Bush construction"
  aus
}
