#' functions for N from Colbourn tables and NISTcat
#'
#' The functions extract the number of runs needed for a specified setting
#' according to the Colbourn tables and according to the NIST catalogues.
#'
#' @rdname CAN
#'
#' @aliases eCAN
#' @aliases N_NISTcat
#' @aliases N_TJcat
#'
#' @usage eCAN(t, k, v)
#' @usage N_NISTcat(t, k, v)
#' @usage N_TJcat(t, k, v)
#'
#' @param t coverage strength
#' @param k number of columns
#' @param v number of levels for each column
#'
#' @returns \code{eCAN} returns a data frame with the smallest known run size
#'  (empirical CAN, based on the Colbourn table)
#'  and the corresponding source entry (Source),
#'  \code{N_NISTcat} returns the run size of the
#'  corresponding catalogued array of the NIST covering array library,
#'  and \code{N_TJcat} the same thing for the small library by
#'  Jorge Torres-Jimenez.
#'  In cases for which there is no entry in the respective
#'  table or library, the returned results are missing values.
#'
#' @examples
#' eCAN(3, 199, 2)
#' N_NISTcat(3, 199, 2)
#' N_TJcat(3, 199, 2) ## equals the best-known array
#'
#' eCAN(4, 199, 2)
#' N_TJcat(4, 199, 2) ## Colbourn table outdated
#'

## obtain known smallest CA sizes from Colbourn catalogue
#' @export
eCAN <- function(t,k,v){
  hilf <- colbournBigFrame[which(colbournBigFrame$t==t & colbournBigFrame$v==v),]
  #  genau <- which(hilf$k==k)
  #  if (length(genau)==1)
  #    return(hilf$N[genau]) else
  hilf <- hilf[which(hilf$k>=k),]
  if (nrow(hilf)==0) return(data.frame(CAN=NA, Source=""))
  zeile <- which.min(hilf$N)
  hilf <- hilf[zeile,c("N","Source")]
  colnames(hilf) <- c("CAN", "Source")
  rownames(hilf) <- NULL
  hilf
}

## obtain size of NIST catalogue entry
#' @export
N_NISTcat <- function(t,k,v){
  hilf <- NISTcat[NISTcat[,"t"]==t & NISTcat[,"k"]==k & NISTcat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(hilf[,"N"])
}

## obtain size of TJ catalogue entry
#' @export
N_TJcat <- function(t,k,v){
  hilf <- TJcat[TJcat[,"t"]==t & TJcat[,"k"]>=k & TJcat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(min(hilf[,"N"]))
}
