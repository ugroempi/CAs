#' Function to juxtapose two strength t CAs
#'
#' Two CAs D1 and D2 with the same strength t, the same number of columns k,
#' and the same numbers of levels in all but at most one column
#' are juxtaposed, such that the levels of one corresponding column
#' in each array are concatenated, resulting in a new strength t CA
#' with the sum of the numbers of levels for that column.
#'
#' @rdname juxtapose
#'
#' @aliases juxtapose
#'
#' @usage juxtapose(D1, D2, jcol=NULL, start0=TRUE, ...)
#'
#' @param D1 a CA with levels as integers 0, ..., v-1 or 1,...,v
#' @param D2 a CA like \code{D1}, same number of columns, same numbers of levels,
#' except perhaps for the \code{jcol}th column; must have the same coding as \code{D3}
#' @param jcol the id of the column to be concatenated; if \code{NULL},
#' it is either the only column id for which the numbers of levels differ, or 1.
#' @param start0 logical: Do integer valued levels start with 0 ? (Otherwise, they must start with 1.)
#' @param ... further arguments to function \code{\link{coverage}}
#'
#' @details
#' Juxtaposing two CAs of a particular strength yields a new CA of the same strength.
#' If \code{D1} and \code{D2} do not have the same strength, the weaker strength prevails.
#'
#' @examples
#' # create a CA(24,3,2^10 4^1) from two CA(12,3,11,2)
#'
#' ## do it by hand
#' A <- paleyHad(11)
#' A2 <- A
#' A2[,11] <- A[,11] + 2
#' coverage(rbind(A,A2),3)
#'
#' ## use the function
#' juxtapose(A, A, jcol=11)
#'
#' ## create a CA(13, 2^10 3^1, 2) by juxtaposing two
#' ##    strength 2 CAs (one with a column that has only 1 level):
#' juxtapose(KSK(11), cbind(KSK(10), 0))
#'    ## note: it is not necessary to specify jcol
#'    ##       because there is only one column with
#'    ##       different number of levels
#'
#' @export
juxtapose <- function(D1, D2, jcol=NULL, start0=TRUE, ...){
  if (!is.matrix(D1)) D1 <- as.matrix(D1)
  if (!is.matrix(D2)) D2 <- as.matrix(D2)
  ll1 <- levels.no(D1)
  ll2 <- levels.no(D2)
  if (is.null(jcol)){
    jcol <- which(!ll1==ll2)
    if (length(jcol)>1)
      stop("The numbers of levels of D1 and D2 must be the same except for at most one column.")
    if (length(jcol)==0) jcol <- 1
  }
  if (!all(ll1[-jcol]==ll2[-jcol]))
    stop("The numbers of levels of D1 and D2 must be the same except for column jcol.")
  if ((!min(D1)==0) && start0) stop("Did you forget to specify start0=FALSE?")
  if ((!start0) && !min(D1)==1) stop("With start0=FALSE, levels must be integers 1,...,v")
  D2[,jcol] <- D2[,jcol] + ll1[jcol]
  aus <- rbind(D1, D2)
  rownames(aus) <- NULL
  aus
}

