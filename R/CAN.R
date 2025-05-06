#' Functions for N and k from Colbourn tables and other libraries
#'
#' The functions extract the number of runs needed for a specified setting
#' according to the Colbourn tables, according to the NIST library of CAs
#' and the Torres-Jimenez library of CAs. Or they extract the number of columns
#' achievable for a specified number of runs.
#'
#' @rdname CAN
#'
#' @aliases eCAN
#' @aliases N_NISTcat
#' @aliases N_TJcat
#' @aliases eCAK
#' @aliases k_NISTcat
#' @aliases k_TJcat
#'
#' @usage eCAN(t, k, v)
#' @usage N_NISTcat(t, k, v)
#' @usage N_TJcat(t, k, v)
#' @usage eCAK(t, N, v)
#' @usage k_NISTcat(t, N, v)
#' @usage k_TJcat(t, N, v)
#'
#' @param t coverage strength
#' @param k number of columns
#' @param v number of levels for each column
#' @param N number of runs
#'
#' @returns \code{eCAN} returns a data frame with the smallest known run size
#'  (empirical CAN, based on the Colbourn table)
#'  and the corresponding source entry (Source),\cr
#'  \code{eCAK} does the same with the largest possible number of columns \code{k},\cr
#'  \code{N_NISTcat} returns the smallest run size for the requested \code{k} of a
#'  catalogued array from the NIST covering array library,\cr
#'  \code{k_NISTcat} returns the maximum column size for the requested \code{N} of a
#'  catalogued array from the NIST covering array library,\cr
#'  and \code{N_TJcat} returns the smallest run size for the requested \code{k} of a
#'  catalogued array from the Jose Torres-Jimenez library (as of Feb 6 2025, small
#'  selection, which is expected to grow),\cr
#'  and \code{k_TJcat} returns the largest number of columns \code{k} for the requested \code{N}
#'  of a catalogued array from the Jose Torres-Jimenez library (as of Feb 6 2025, small
#'  selection, which is expected to grow)
#'  .
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
#' eCAK(2,11,3)  ## also obtained by CAEX
#' k_NISTcat(2,11,3)
#' k_TJcat(2,11,3) ## equals the best-known array
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

#' @export
eCAK <- function(t,N,v){
  hilf <- colbournBigFrame[which(colbournBigFrame$t==t & colbournBigFrame$v==v),]
  #  genau <- which(hilf$k==k)
  #  if (length(genau)==1)
  #    return(hilf$N[genau]) else
  hilf <- hilf[which(hilf$N<=N),]
  if (nrow(hilf)==0) return(data.frame(CAK=NA, Source=""))
  zeile <- which.max(hilf$k)
  hilf <- hilf[zeile,c("k","Source")]
  colnames(hilf) <- c("CAK", "Source")
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

#' @export
k_NISTcat <- function(t,N,v){
  hilf <- NISTcat[NISTcat[,"t"]==t & NISTcat[,"N"]<=N & NISTcat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(max(hilf[,"k"]))
}

## obtain size of TJ catalogue entry
#' @export
N_TJcat <- function(t,k,v){
  hilf <- TJcat[TJcat[,"t"]==t & TJcat[,"k"]>=k & TJcat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(min(hilf[,"N"]))
}

#' @export
k_TJcat <- function(t,N,v){
  hilf <- TJcat[TJcat[,"t"]==t & TJcat[,"N"]<=N & TJcat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(max(hilf[,"k"]))
}
