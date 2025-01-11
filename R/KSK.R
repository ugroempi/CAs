#' function for optimal 2-level CAs of strength 2
#'
#' implements the construction by Kleitman and Spencer (1973) and Katona (1973),
#' which have optimal size
#'
#' @rdname KSK
#' @aliases k_KSK
#' @aliases N_KSK
#'
#' @usage KSK(k=NULL, N=NULL)
#' @usage k_KSK(N)
#' @usage N_KSK(k)
#'
#'
#' @param k number of columns
#' @param N number of runs
#'
#' @return \code{KSK} returns a matrix with \code{k} columns and \code{N} runs and levels 0 and 1
#' in each column, with coverage strength 2.
#'
#' \code{N_KSK} and \code{k_KSK} return the minimum \code{N} for a given \code{k}
#' or the maximum \code{k} for a given \code{N}, respectively.
#'
#' @examples
#' ## smallest design for 10 columns
#' KSK(k=10)
#' ## largest design for 6 rows
#' ca6.2.10.2 <- KSK(N=6)
#' coverage(ca6.2.10.2, 2)
#' coverage(ca6.2.10.2, 3)
#' ## largest k for 10 rows
#' k_KSK(10)
#' ## smallest N for 2800 columns
#' N_KSK(2800)
#' ## design for 8 columns in 6 rows (not necessarily best subset of columns)
#' KSK(8,6)
#' ## design for 2 columns (full factorial needed)
#' KSK(2)
#'

N_KSK <- function(k) KSKnk$N[min(which(KSKnk$maxk>=k))]

k_KSK <- function(N) choose(N-1, ceiling(N/2))

KSK <- function(k=NULL, N=NULL){
  ## k is the number of columns
  ## N is the number of rows
  ## if one is NULL, it is calculated from the other
  ## issue: the design for k < maxk is not necessarily the best choice of columns
  ## KSKnk is a data frame in sysdata.rda and holds the maximum number of columns
  ##    that can be done in a given number of runs.
  if (is.null(k) && is.null(N)) stop("At least one of k and N must be specified")
  if (is.null(k)) k <- k_KSK(N)
  if (is.null(N)) N <- N_KSK(k)
  pos1 <- nchoosek(N-1, ceiling(N/2))
  CA <- matrix(0, N, k)
  for (j in 1:k){
    CA[pos1[,j]+1,j] <- 1
  }
  CA
}

