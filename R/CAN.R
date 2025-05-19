#' Functions for N and k from Colbourn tables, other libraries, and implemented constructions
#'
#' The functions extract the number of runs needed for a specified setting
#' according to the Colbourn tables, according to the NIST library of CAs
#' and the Torres-Jimenez library of CAs, as well as according to implemented constructions
#' and available catalogues. Or they extract the number of columns
#' achievable for a specified number of runs.
#'
#' @rdname CAN
#'
#' @aliases Ns
#' @aliases eCAN
#' @aliases N_NISTcat
#' @aliases N_TJcat
#' @aliases N_WKScat
#' @aliases N_CKRScat
#' @aliases N_DWYERcat
#' @aliases N_PALEYcat
#' @aliases N_CYCLOTOMYcat
#' @aliases N_CAEX
#' @aliases ks
#' @aliases eCAK
#' @aliases k_NISTcat
#' @aliases k_TJcat
#' @aliases k_WKScat
#' @aliases k_CKRScat
#' @aliases k_DWYERcat
#' @aliases k_PALEYcat
#' @aliases k_CYCLOTOMYcat
#' @aliases k_CAEX
#'
#' @usage Ns(t, k, v)
#' @usage eCAN(t, k, v)
#' @usage N_NISTcat(t, k, v)
#' @usage N_TJcat(t, k, v)
#' @usage N_WKScat(t=6, k, v=2)
#' @usage N_CKRScat(t, k, v)
#' @usage N_DWYERcat(t, k, v)
#' @usage N_PALEYcat(t, k, v=2)
#' @usage N_CYCLOTOMYcat(t, k, v)
#' @usage N_CAEX(t=2, k, v=3)
#' @usage ks(t, N, v)
#' @usage eCAK(t, N, v)
#' @usage k_NISTcat(t, N, v)
#' @usage k_TJcat(t, N, v)
#' @usage k_WKScat(t=6, N, v=2)
#' @usage k_CKRScat(t, N, v)
#' @usage k_DWYERcat(t, N, v)
#' @usage k_PALEYcat(t, N, v=2)
#' @usage k_CYCLOTOMYcat(t, N, v)
#' @usage k_CAEX(t=2, N, v=3)
#'
#' @param t coverage strength
#' @param k number of columns
#' @param v number of levels for each column
#' @param N number of runs
#'
#' @details
#' Functions \code{Ns} and \code{ks} take into account all available
#' catalogues and constructions, for which sizes can be easily provided; their
#' scope will grow with time. Specific functions for run sizes \code{N} and
#' numbers of columns \code{k} can also be inspected separately.
#' The functions typically look for a value of \code{k} for \emph{at most}
#' \code{N} runs, but a value of \code{N} for exactly \code{k} columns.
#' Most functions look for
#' exact strength \code{t}, except for those for \code{WKScat}, which also return a
#' value for at most strength \code{t}, with a message.
#'
#' At present, the individual functions,
#' except for \code{eCAN} and \code{eCAK}, return a number only.
#' This may change in the future.
#'
#' For the \code{TJcat}-related functions, a message provides the command for
#' creating the CA, if a creation from within the package is possible (with one of
#' functions \code{CAEX} or \code{SCA_Bose}).
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
#' N_CYCLOTOMYcat(3, 45, 3) ## t=3 and v=3 not implemented
#' N_CYCLOTOMYcat(4, 45, 3) ## quite large!
#' k_CYCLOTOMYcat(4, 669, 3) ## can accommodate 224 columns
#' N_CYCLOTOMYcat(4, 5, 3) ## is the smallest *implemented*
#'                         ## cyclotomy strength 4 CA for 3 levels
#' N_CYCLOTOMYcat(4, 670, 3) ## the next smallest is for 1051 columns
#' eCAN(4, 670, 3)           ## which is a current optimum
#'
#' # overview of available constructions
#' Ns(2, 45, 2)  ## function KSK is known to yield the overall optimum
#' Ns(3, 45, 2)  ## the best array is in the DWYER-catalogue
#'    # https://github.com/aadwyer/CA_Database/blob/main/Repository/CA/CA_26_3_46_2.txt
#' Ns(3, 400, 2) ## at present, the best easily available array
#'               ## is in the NIST catalogue
#'    # https://math.nist.gov/coveringarrays/ipof/cas/t=3/v=2/ca.3.2%5E400.txt.zip
#' Ns(2, 800, 3) ## function CAEX produces the best available array
#' Ns(3, 200, 3) ## the respective construction has not yet been implemented
#'    # ## at present, the best easily available array
#'               ## is in the NIST catalogue
#'    # https://math.nist.gov/coveringarrays/ipof/cas/t=3/v=3/ca.3.3%5E200.txt.zip
#' Ns(4, 150, 3) ## cyclotomy is very competitive
#' Ns(4, 150, 2) ## the best easily available CA is NIST, PALEY is not much worse
#' Ns(4, 50, 2) ## Paley and Cyclotomy are the same, and both optimal
#' Ns(4, 50, 3) ## the best easily available CA is
#'    # https://github.com/aadwyer/CA_Database/blob/main/Repository/CA/CA_507_4_53_3.txt
#' Ns(5, 50, 2) ## Paley is optimal
#' Ns(5, 500, 2) ## at present, Cyclotomy is competitive - best design will be implemented
#' Ns(5, 50, 3) ## the best easily available CA is
#'    # https://github.com/aadwyer/CA_Database/blob/main/Repository/CA/CA_2067_5_50_3.txt
#' Ns(5, 500, 3) ## Cyclotomy is optimal
#' Ns(6, 50, 2) ## WKS is optimal
#'    # WKS_CAs['50']
#' Ns(6, 200, 2) ## Paley is optimal
#' Ns(6, 2000, 2) ## Cyclotomy is optimal
#' Ns(6, 200, 3) ## not easily available at present
#'

## provide sizes of implemented methods and unimplemented catalogues
#' @export
Ns <- function(t, k, v){
  aus <- c(
    PALEY=N_PALEYcat(t,k,v),
    CAEX=N_CAEX(t,k,v),
    CYCLOTOMY=N_CYCLOTOMYcat(t,k,v),
    CKRS=N_CKRScat(t,k,v),
    WKS=N_WKScat(t,k,v),
    CS_MS=N_CS_MS(k,v),
    CS_LCDST=N_CS_LCDST(k,v),
    DWYER=N_DWYERcat(t,k,v),
    NIST=N_NISTcat(t,k,v),
    TJ=N_TJcat(t,k,v),
    eCAN=eCAN(t,k,v)[[1]])
  aus <- aus[!is.na(aus)]
  if (t==2 && v==2) aus <- c(KSK=N_KSK(k), aus)
  if (length(aus)==0){
    message("no solution found")
    return(numeric(0))
  }
  aus
}

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

N_CAEX <- function(t=2,k,v=3){
  if (t>2 || !v==3) return(NA)
  N_TJcat(t,k,v)
}

## obtain size of WKS catalogue entry
#' @export
N_WKScat <- function(t=6,k,v=2){
  if (!v==2) return(NA)
  if (t>6 || k>max(WKScat$k)) return(NA)
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

## obtain size of CYCLOTOMY catalogue entry
#' @export
N_CYCLOTOMYcat <- function(t,k,v){
  hilf <- CYCLOTOMYcat[CYCLOTOMYcat[,"t"]==t & CYCLOTOMYcat[,"k"]>=k &
                         CYCLOTOMYcat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(min(hilf[,"N"]))
}

## obtain size of PAYLEY catalogue entry
#' @export
N_PALEYcat <- function(t,k,v=2){
  hilf <- PALEYcat[PALEYcat[,"t"]==t & PALEYcat[,"k"]>=k & PALEYcat[,"v"]==v,,drop=FALSE]
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

ks <- function(t, N, v){
  aus <- c(
    PALEY=k_PALEYcat(t,N,v),
    CAEX=k_CAEX(t,N,v),
    CYCLOTOMY=k_CYCLOTOMYcat(t,N,v),
    CKRS=k_CKRScat(t,N,v),
    WKS=k_WKScat(t,N,v),
    CS_MS=k_CS_MS(N,v),
    CS_LCDST=k_CS_LCDST(N,v),
    DWYER=k_DWYERcat(t,N,v),
    NIST=k_NISTcat(t,N,v),
    TJ=k_TJcat(t,N,v),
    eCAK=eCAK(t,N,v)[[1]])
  aus <- aus[!is.na(aus)]
  if (t==2 && v==2) aus <- c(KSK=k_KSK(N), aus)
  if (length(aus)==0){
    message("no solution found")
    return(numeric(0))
  }
  aus
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

k_CAEX <- function(t=2,N,v=3){
  if (t>2 || !v==3) return(NA)
  k_TJcat(t,N,v)
}

#' @export
k_DWYERcat <- function(t,N,v){
  hilf <- DWYERcat[DWYERcat[,"t"]==t & DWYERcat[,"N"]<=N & DWYERcat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(max(hilf[,"k"]))
}

#' @export
k_CYCLOTOMYcat <- function(t,N,v){
  hilf <- CYCLOTOMYcat[CYCLOTOMYcat[,"t"]==t & CYCLOTOMYcat[,"N"]<=N & CYCLOTOMYcat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(max(hilf[,"k"]))
}

#' @export
k_PALEYcat <- function(t,N,v=2){
  hilf <- PALEYcat[PALEYcat[,"t"]==t & PALEYcat[,"N"]<=N & PALEYcat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(max(hilf[,"k"]))
}

#' @export
k_WKScat <- function(t=6,N,v=2){
  if (!v==2) return(NA)
  if (t>6) return(NA)
  if (t<6) message("The WKS catalogue has strength 6 CAs only.")
  hilf <- WKScat[WKScat[,"N"]<=N,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(max(hilf[,"k"]))
}

#' @export
k_CKRScat <- function(t,N,v){
  hilf <- CKRScat[CKRScat[,"t"]==t & CKRScat[,"N"]<=N & CKRScat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(max(hilf[,"k"]))
}

