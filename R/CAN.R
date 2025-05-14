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
#' @aliases N_WKScat
#' @aliases N_CKRScat
#' @aliases N_DWYERcat
#' @aliases eCAK
#' @aliases k_NISTcat
#' @aliases k_TJcat
#' @aliases k_WKScat
#' @aliases k_CKRScat
#' @aliases k_DWYERcat
#'
#' @usage eCAN(t, k, v)
#' @usage N_NISTcat(t, k, v)
#' @usage N_TJcat(t, k, v)
#' @usage N_WKScat(t, k, v)
#' @usage N_CKRScat(t, k, v)
#' @usage N_DWYERcat(t, k, v)
#' @usage eCAK(t, N, v)
#' @usage k_NISTcat(t, N, v)
#' @usage k_TJcat(t, N, v)
#' @usage k_WKScat(t, N, v)
#' @usage k_CKRScat(t, N, v)
#' @usage k_DWYERcat(t, N, v)
#'
#' @param t coverage strength
#' @param k number of columns
#' @param v number of levels for each column
#' @param N number of runs
#'
#' @details
#' The functions typically look for a value of \code{k} for at most \code{N} runs,
#' but a value of \code{N} for exactly \code{k} columns. Most functions look for
#' exact strength \code{t}, except for those for \code{WKScat}, which also return a
#' value for at most strength \code{t}, with a message.
#'
#' At present, all functions
#' except for \code{eCAN} and \code{eCAK} return a number only.
#' This may change in the future.
#'
#' For the \code{TJcat}-related functions, a message provides the command for
#' creating the CA, if a creation from within the package is possible (with one of
#' functions \code{CAEX} or \code{SCA_Bose}).
#'
#' The catalogue objects are internal;
#' they can thus be accessed by using \code{CAs:::TJcat} etc.
#' This may change in the future.
#'
#'
#' @returns \code{eCAN} returns a data frame with the smallest known run size
#'  (empirical CAN, based on the Colbourn table)
#'  and the corresponding source entry (Source),\cr
#'  \code{eCAK} does the same with the largest possible number of columns \code{k},\cr
#'  \code{N_NISTcat} returns the smallest run size for the requested \code{k} of a
#'  catalogued array from the publicly available NIST covering array library,\cr
#'  \code{k_NISTcat} returns the maximum column size for the requested \code{N} of a
#'  catalogued array from the publicly available NIST covering array library,\cr
#'  and \code{N_TJcat} returns the smallest run size for the requested \code{k} of a
#'  catalogued array from the Jose Torres-Jimenez library (as of Feb 6 2025, small
#'  selection, which is expected to grow; partly available in this package),\cr
#'  and \code{k_TJcat} returns the largest number of columns \code{k} for the requested \code{N}
#'  of a catalogued array from the Jose Torres-Jimenez library (as of Feb 6 2025, small
#'  selection, which is expected to grow; partly available in this package)\cr
#'  The analogous functionality holds for \code{N_} and \code{k_} functions for
#'  \code{DWYERcat}, \code{CKRScat} and \code{WKScat}.
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
  else {
    if (t==2 && v==3){
      if (k==4)
        message("The CA can be produced with command SCA_Bose(3).")
      else
        message("The CA can be produced with command CAEX(",k,").")
    }
    return(min(hilf[,"N"]))
  }
}

## obtain size of WKS catalogue entry
#' @export
N_WKScat <- function(t=6,k,v=2){
  if (!v==2) return(NA)
  if (t>6) return(NA)
  if (t<6) message("The WKS catalogue has strength 6 CAs only.")
  hilf <- WKScat[WKScat[,"k"]>=k,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(min(hilf[,"N"]))
}

## obtain size of DWYER catalogue entry
#' @export
N_DWYERcat <- function(t,k,v){
  hilf <- DWYERcat[DWYERcat[,"t"]==t & DWYERcat[,"k"]>=k & DWYERcat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(min(hilf[,"N"]))
}

## obtain size of CKRS catalogue entry
#' @export
N_CKRScat <- function(t,k,v){
  hilf <- CKRScat[CKRScat[,"t"]==t & CKRScat[,"k"]>=k & CKRScat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(min(hilf[,"N"]))
}

#' @export
k_TJcat <- function(t,N,v){
  hilf <- TJcat[TJcat[,"t"]==t & TJcat[,"N"]<=N & TJcat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else{
    k <- max(hilf[,"k"])
    if (t==2 && v==3){
      if (N < 11)
        message("The CA can be produced with command SCA_Bose(3).")
      else
        message("The CA can be produced with command CAEX(",k,").")
    }
  return(k)
  }
}

#' @export
k_DWYERcat <- function(t,N,v){
  hilf <- DWYERcat[DWYERcat[,"t"]==t & DWYERcat[,"N"]<=N & DWYERcat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(max(hilf[,"k"]))
}

#' @export
k_WKScat <- function(t,N,v){
  if (!v==2) return(NA)
  if (t>6) return(NA)
  if (t<6) message("The WKS catalogue has strength 6 CAs only.")
  hilf <- WKScat[WKScat[,"N"]<=N,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(max(hilf[,"k"]))
}

k_CKRScat <- function(t,N,v){
  hilf <- CKRScat[CKRScat[,"t"]==t & CKRScat[,"N"]<=N & CKRScat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(max(hilf[,"k"]))
}

