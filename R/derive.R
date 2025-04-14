#' Function to derive a strength t+1 CA into a smaller strength t array
#'
#' A CA D of strength t+1 is reduced to a single level of a particular column,
#' which is subsequently omitted from the array.
#' The resulting CA has strength t.
#'
#' @rdname derive
#'
#' @aliases derive
#'
#' @usage derive(D, dcol=1, dval=NULL, start0=TRUE, ...)
#'
#' @param D a CA with levels as integers 0, ..., v-1 or 1,...,v
#' @param dcol the id of the column to be deleted (default: 1)
#' @param dval the value to be used for filtering (default: the first value that has minimum frequency in column \code{jcol})
#' @param start0 logical: Do integer valued levels start with 0 ? (Otherwise, they must start with 1.)
#' @param ... currently unused
#'
#' @details
#' The name \code{derive} was taken from Zhang et al. (2014).
#'
#'
#' @examples
#' # create a CA(6,2,10,2) from a CA(12,3,11,2)
#'
#' A <- paleyHad(11)  ## a CA(12,3,11,2)
#'
#' sixdefault <- derive(A)  ## first column of A was 0,
#'                          ## first column omitted
#' dim(sixdefault)
#' coverage(sixdefault, 2)
#'
#' ## select all the runs with 11th column at level 1
#' ## keep first 10 columns for these
#' six <- derive(A, dcol=11, dval=1)
#' dim(six)
#' ## strength 3 was reduced to strength 2
#' coverage(six, 2)
#' eCAN(2,10,2)  ## the array is optimal
#'
#'
#' @export
derive <- function(D, dcol=1, dval=NULL, start0=TRUE, ...){
  if (!is.matrix(D)) D <- as.matrix(D)
  ll <- levels.no(D)
  if (!dcol %in% 1:ncol(D)) stop("dcol must be a column number for D")
  if (is.null(dval)){
    tab <- table(D[,dcol])
    dval <- as.numeric(names(tab)[which.min(tab)])
  }
  stopifnot(dval %in% unique(D[,dcol]))
  if ((!min(D)==0) && start0) stop("Did you forget to specify start0=FALSE?")
  if ((!start0) && !min(D)==1) stop("With start0=FALSE, levels must be integers 1,...,v")
  aus <- D[D[,dcol]==dval,][,-dcol]
  dimnames(aus) <- NULL
  aus
}

