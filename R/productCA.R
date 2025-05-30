#' Direct product and generalized direct product for two CAs
#'
#' Functions to construct a CA from a (generalized) direct product of two ingoing CAs
#'
#' @rdname productCA
#'
#' @aliases productCA
#' @aliases productCA_raw
#'
#' @usage productCA(D1, D2, check=TRUE, dupremove=TRUE, generalized=TRUE, ...)
#' @usage productCA_raw(D1, D2)
#'
#' @param D1 an N x k CA of strength 2 with v levels
#' @param D2 an M x l CA of strength 2 with v levels
#' @param check logical; if TRUE, checks required strength of ingoing CAs (may substantially increase run time for large CAs)
#' @param dupremove logical; if TRUE, removes duplicated rows
#' @param generalized logical; if TRUE, exploits constant rows, which implies that function \code{\link{maxconstant}} is used to maximize constant rows,
#' which can reduce \code{N+M} by up to \code{v} (and at least by two); the price is run time.
#' @param ... further arguments to function \code{\link{maxconstant}}, e.g., \code{one_is_enough}, which is per default FALSE
#'
#' @section Details:
#' Function \code{productCA} yields a CA(N + M -  min(v, s+r), 2, kl, v) from a product
#' construction,\cr
#' for which \code{D1} can be rearranged to have s constant rows \cr
#' and \code{D2} can be rearranged to have r constant rows \cr
#' as proposed in Corollary 4.3 of Colbourn and Torres-Jimènez (2013)\cr
#' (removal of duplicates may reduce run size further, for very poor ingoing arrays).\cr
#' For large and efficient ingoing arrays, it may be wise to specify \code{one_is_enough=TRUE}
#' (as part of \code{...});
#' in that case, the search for the maximum number of constant rows is skipped,
#' and the resulting run size is N + M - 2.
#'
#' Function \code{productCA} also serves as work horse function for combining two CAs in some instances of
#' the \code{\link{CAEX}} construction.
#'
#' Function \code{productCA_raw} is mainly meant for calculation of interim results in other
#' functions, e.g., \code{\link{productPCA}}.
#' It works on any pair of numeric matrices that do not need to be integer-valued or CAs (see Examples Section).
#'
#' @returns Function \code{productCA} returns a CA(N + M -  min(v, s+r), 2, kl, v) (matrix of class \code{ca}).\cr
#' Function \code{productCA_raw} returns the simple (N + M) x (k*l) direct product of its two arguments.
#'
#' @references Colbourn and Torres-Jimenez (2013)
#'
#' @examples
#' #########################################
#' ## productCA
#' #########################################
#' ###  (N+M-2 x k*l matrix)
#' A <- cyc(19,2)
#' B <- KSK(k=12)
#' dim(A); dim(B)
#' E <- productCA(A, B)
#' dim(E)
#' coverage(E, 2)
#' ## 19 + 7 - 2 runs, twice as large as the optimal 12
#' eCAN(2,228,2)
#'
#' A <- lhs::createBose(4, 5, bRandom=FALSE)
#' dim(A)
#' E <- productCA(A, A)
#' coverage(E, 2)
#' dim(E)  ## 30 runs with 25 columns, close to optimal
#' eCAN(2, 25, 4) ## 29 runs
#' eCAK(2, 30, 4) ## 30 columns
#'
#' ## from two CAs with many constant rows
#' ## (not useful for this particular A, for demo purposes only)
#' A <- lhs::createBush(4, 5, bRandom=FALSE)
#' dim(productCA(A,A)) ## using default behavior
#'                      ## that exploits maximum constant rows
#'                      ## reduces N+M by up to v runs
#' dim(productCA(A,A, one_is_enough=TRUE))
#'                      ## preventing function maxconstant from conducting a
#'                      ## clique search, reduces N+M by 2 runs
#'
#' ## the function also works for designs with flexible values (which are denoted by NA)
#' ## productCA
#' D <- productCA(CAEX(N=16), CAEX(N=16))
#' dim(D)
#' eCAN(2,441,3) ## slightly worse than the best CA, which is obtained as CAEX(N=28)
#' eCAK(2,30,3) ## a lot more columns are possible in 30 runs
#'
#' ### productCA_raw works on arbitrary numeric matrices
#' ### and applies the simple direct product
#' A <- matrix(1:12, nrow=4)
#' B <- matrix(10*(1:6), nrow=3)
#' productCA_raw(A, B)
#'

#' @export
productCA_raw <- function(D1, D2){
  ## function for use in productPCA
  ## skips all the checks, except for matrix
  ## omit all constant row processing
  stopifnot(is.matrix(D1), is.matrix(D2))
  stopifnot(is.numeric(D1), is.numeric(D2))
  rbind(
        do.call(cbind, rep(list(D1), ncol(D2))),
        kronecker(D2, matrix(1, 1, ncol(D1)))
        )
}

#' @export
productCA <- function(D1, D2, check=TRUE, dupremove=TRUE, generalized=TRUE, ...){
  Call <- sys.call()
  stopifnot(is.matrix(D1)); stopifnot(is.matrix(D2))
  v <- levs1 <- unique(levels.no.NA(D1))
  levs2 <- unique(levels.no.NA(D2))
  ## at present, only for uniform CAs
  stopifnot(length(levs1)==1);stopifnot(length(levs2)==1)
  ## both CAs need to have the same number of columns
  stopifnot(levs1==levs2)
  # infer start0
  start0 <- (min(D1, na.rm=TRUE)==0)
  if(check){
    stopifnot(all(coverage(D1, 2, start0=start0)==1))
    stopifnot(all(coverage(D2, 2, start0=start0)==1))
  }
  k <- ncol(D1);
  l <- ncol(D2)
  if (generalized){
    D1 <- maxconstant(D1, verbose=2, ...);
    nc1 <- length(attr(D1, "constant_rows")$row_set_list)
    if (nc1 == v){
      D1 <- D1[-(1:v),,drop=FALSE]
      nc2 <- sum(apply(D2, 1, function(obj) length(unique(obj))==1))
    }else{ ## nc1 < v
      D2 <- maxconstant(D2, verbose=2, ...)
      nc2 <- length(attr(D2, "constant_rows")$row_set_list)

      ## ensure disjoint constant rows
      constvals1 <- D1[1:nc1, 1, drop=TRUE]
      constvals2 <- D2[1:nc2, 1, drop=TRUE]
      targetconstvals2 <- setdiff((1:v) - as.numeric(start0), constvals1)
      if (length(union(constvals1, constvals2))==min(v, nc1+nc2)){
        ## maximum possible disjoint rows already
        D1 <- D1[-(1:nc1),]
        D2 <- D2[-(which(constvals2 %in% targetconstvals2)),]
      }else{
        D1 <- D1[-(1:nc1),]   ## assuming these rows are disjoint
        ## achieve disjoint constant rows in D2
        for (r in 1:min(nc2,length(targetconstvals2)))
          D2 <- swapvals(D2, 1:l, D2[r,1], targetconstvals2[r])
        D2 <- D2[-(1:min(nc2,length(targetconstvals2))),]
      }
    }## end of nc1 < v
    }## end of if (generalized)

  ## product with up to v disjoint constant rows removed
  aus <- productCA_raw(D1, D2)
  if (dupremove) aus <- unique(aus)
  aus <- CA_to_PCA(aus)
  attr(aus, "Call") <- Call
  aus
}
