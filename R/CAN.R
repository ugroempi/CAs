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
#' @aliases Ns_derive
#' @aliases Ns_fuse
#' @aliases Ns_CK_doubling
#' @aliases Ns_productCA
#' @aliases eCAN
#' @aliases N_NISTcat
#' @aliases N_TJcat
#' @aliases N_WKScat
#' @aliases N_CKRScat
#' @aliases N_DWYERcat
#' @aliases N_PALEYcat
#' @aliases N_CYCLOTOMYcat
#' @aliases N_CAEX
#' @aliases N_recBoseCA
#' @aliases ks
#' @aliases ks_productCA
#' @aliases eCAK
#' @aliases k_NISTcat
#' @aliases k_TJcat
#' @aliases k_WKScat
#' @aliases k_CKRScat
#' @aliases k_DWYERcat
#' @aliases k_PALEYcat
#' @aliases k_CYCLOTOMYcat
#' @aliases k_CAEX
#' @aliases k_recBoseCA
#'
#' @usage Ns(t, k, v)
#' @usage Ns_derive(t, k, v)
#' @usage Ns_fuse(t, k, v, maxfuse = 1)
#' @usage Ns_CK_doubling(t=3, k, v)
#' @usage Ns_productCA(t=2, k, v)
#' @usage eCAN(t, k, v)
#' @usage N_NISTcat(t, k, v)
#' @usage N_TJcat(t, k, v)
#' @usage N_WKScat(t=6, k, v=2)
#' @usage N_CKRScat(t, k, v)
#' @usage N_DWYERcat(t, k, v)
#' @usage N_PALEYcat(t, k, v=2)
#' @usage N_CYCLOTOMYcat(t, k, v)
#' @usage N_CAEX(t=2, k, v=3)
#' @usage N_recBoseCA(t=2, k, v, type="PCA")
#' @usage ks(t, N, v)
#' @usage Ns(t, k, v)
#' @usage ks_productCA(t=2, N, v)
#' @usage eCAK(t, N, v)
#' @usage k_NISTcat(t, N, v)
#' @usage k_TJcat(t, N, v)
#' @usage k_WKScat(t=6, N, v=2)
#' @usage k_CKRScat(t, N, v)
#' @usage k_DWYERcat(t, N, v)
#' @usage k_PALEYcat(t, N, v=2)
#' @usage k_CYCLOTOMYcat(t, N, v)
#' @usage k_CAEX(t=2, N, v=3)
#' @usage k_recBoseCA(t=2, N, v, type="PCA")
#'
#' @param t coverage strength
#' @param k number of columns
#' @param v number of levels for each column
#' @param maxfuse integer number of levels to fuse (i.e., difference between number of levels before and after fusing)
#' @param N number of runs
#' @param type character string: \code{"PCA"} or \code{"CA"}
#'
#' @details
#' Functions \code{Ns} and \code{ks} take into account all available
#' catalogues and constructions, for which sizes can be easily provided; their
#' scope will grow with time. They provide the respective run sizes for
#' a requested number of columns \code{k} (\code{Ns}) or the affordable numbers
#' of columns from an affordable run size \code{N} (\code{ks}).\cr
#' Functions \code{Ns_derive}, \code{Ns_fuse}, \code{Ns_CK_doubling}, \code{Ns_productCA}
#' provide run sizes from deriving, fusing, doubling or taking the product with
#' the smallest possible strength 2 CA for the requested \code{v}
#' in comparison to run sizes without applying such techniques. In case of doubling,
#' the function is restricted to strength $t=3$ (the only strength for which
#' Chateauneuf-Kreher doubling works) with $v=2,3$, for which a natural choice of a
#' the strength 2 design \code{D2} for \code{\link{CK_doubling}} is known (implemented
#' in functions \code{\link{KSK}} and \code{\link{CAEX}}, respectively).\cr
#' Functions \code{ks_derive} and \code{ks_fuse} provide affordable
#' numbers of columns from an affordable run size \code{N};
#' for doubling, this is more complicated and has not yet been implemented.
#'
#' Specific functions for run sizes \code{N} and
#' numbers of columns \code{k} can also be inspected separately.\cr
#' The functions typically look for a value of \code{k} for \emph{at most}
#' \code{N} runs, and a value of \code{N} for at least \code{k} columns.\cr
#' The functions look for the exact strength \code{t}; this implies that a CA for a higher strength may exist,
#' where none is implemented for a lower strength.
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
#'  \code{DWYERcat}, \code{CKRScat}, \code{WKScat} and \code{CYCLOTOMYcat}.
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
#'    # readable with internet connection
#'    # pathGH <- "https://raw.githubusercontent.com/aadwyer/CA_Database/main/Repository/CA"
#'    # the instruction line starts with a comment character and is ignored, therefore ninstruct=0
#'    # D <- readCA(paste0(pathGH, "/CA_26_3_46_2.txt"), ignore.chars=c("]", "["),
#'    #        ninstruct=0, sep=",")
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
#' Ns_CK_doubling(3, 22, 2)
#' ## with doubling a Paley design for ceiling(22/2) columns,
#' ## the optimum run size can be achieved
#'
#' Ns_CK_doubling(3, 12, 3)
#' ## with doubling a Paley design for ceiling(22/2) columns,
#' ## the optimum run size can be achieved

## provide sizes of implemented methods and unimplemented catalogues
#' @export
Ns <- function(t, k, v){
  suppressMessages(
    aus <- c(
    PALEY=N_PALEYcat(t,k,v),
    CAEX=N_CAEX(t,k,v),
    CYCLOTOMY=N_CYCLOTOMYcat(t,k,v),
    CKRS=N_CKRScat(t,k,v),
    recBoseCA_PCA=unname(N_recBoseCA(t,k,v,type="PCA")),
    recBoseCA_CA=unname(N_recBoseCA(t,k,v,type="CA")),
    WKS=N_WKScat(t,k,v),
    CS_MS=N_CS_MS(k,v),
    CS_LCDST=N_CS_LCDST(k,v),
    DWYER=N_DWYERcat(t,k,v),
    NIST=N_NISTcat(t,k,v),
    TJ=N_TJcat(t,k,v),
    eCAN=eCAN(t,k,v)[[1]])
    )
  aus <- aus[!is.na(aus)]
  if (t==2 && v==2) aus <- c(KSK=N_KSK(k), aus)
  if (length(aus)==0){
    Nderive <- Ns_derive(t, k, v)
    if (!is.na(Nderive)){
      message("no solution found, but deriving yields a solution:")
      aus <- Nderive["N",]
      attr(aus, "detail") <- Nderive
      return(aus)
    }
    message("no solution found")
    return(numeric(0))
  }
  aus
}

## obtain possibilities for derivation sources
#' @export
Ns_derive <- function(t, k, v){
  Nswoderive <- Ns(t, k, v)
  if (t>=6){
    message("derive could not be used for strength 6")
    return(Nswoderive)
  }
  fromNs <- Ns(t+1, k+1, v)
  Ns <- sapply(fromNs, function(obj) deriveto(obj, k+1, v)[1])
  Ns <- as.data.frame(rbind('t-from'=t+1, 'k-from'=k+1, 'N-from'=fromNs, t=t, 'N'=Ns, k=k))
  aus <- dplyr::bind_rows(Nswoderive, Ns, .id="label")
  if (any(aus[1,] > aus[6,], na.rm=TRUE)) return(aus) else {
    message("no derived design was better than direct constructions")
    return(Nswoderive)
  }
}

## obtain possibilities for derivation sources
#' @export
Ns_fuse <- function(t, k, v, maxfuse=1){
  (hilf <- lapply(1:maxfuse,
                 function(obj){
                   hilf <- Ns(t, k, v+obj)
                   aus <- hilf - 2*obj
                   attr(aus, "detail") <- rbind('v-from'=v+obj, 'N_from'=hilf)
                   aus
                 }))
}

## obtain sizes for CK doubling for v=2 and v=3
#' @export
Ns_CK_doubling <- function(t=3, k, v){
  if (!v %in% c(2,3)) return(NA)
  if(!t==3) return(NA)
  khilf <- ceiling(k/2)
  hilf <- Ns(3, khilf, v)
  if (v==3) plus <- suppressMessages(N_CAEX(k=khilf)) else
    plus <- N_KSK(khilf)
  doubled <- as.data.frame(rbind(hilf + plus))
  direct <- as.data.frame(rbind(Ns(t, k, v)))
  dplyr::bind_rows(direct=direct, doubled=doubled, .id="way")
}

## obtain sizes for product construction
## with large D1 and small D2
#' @export
Ns_productCA <- function(t=2, k, v){
  if (!t==2) return(NA)
  M <- v^2
  if (v %in% primedat$q)
    l <- v+1
  else
    l <- 3 ## latin square
  khilf <- ceiling(k/l)
  hilf <- Ns(2, khilf, v)
  productCA <- as.data.frame(rbind(hilf + M - 2))
  suppressMessages(direct <- as.data.frame(rbind(Ns(t, k, v))))
  dplyr::bind_rows(direct=direct, productCA=productCA, .id="way")
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

#' @export
N_CAEX <- function(t=2,k,v=3){
  if (t>2 || !v==3) return(NA)
  N_TJcat(t,k,v)
}

#' @export
N_recBoseCA <- function(t=2,k,v, type="PCA"){
  if (t>2 || !v%in% primedat$q) return(NA)
  suppressMessages(N_d_recursiveBose(v, k, type=type))[1]
}

## obtain size of WKS catalogue entry
#' @export
N_WKScat <- function(t=6,k,v=2){
  if (!v==2) return(NA)
  if (!t==6 || k>max(WKScat$k)) return(NA)
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

#' @export
ks <- function(t, N, v){
  suppressMessages(
    aus <- c(
    PALEY=k_PALEYcat(t,N,v),
    CAEX=k_CAEX(t,N,v),
    CYCLOTOMY=k_CYCLOTOMYcat(t,N,v),
    CKRS=k_CKRScat(t,N,v),
    recBoseCA_PCA=unname(k_recBoseCA(t,N,v,type="PCA")),
    recBoseCA_CA=unname(k_recBoseCA(t,N,v,type="CA")),
    WKS=k_WKScat(t,N,v),
    CS_MS=k_CS_MS(N,v),
    CS_LCDST=k_CS_LCDST(N,v),
    DWYER=k_DWYERcat(t,N,v),
    NIST=k_NISTcat(t,N,v),
    TJ=k_TJcat(t,N,v),
    eCAK=eCAK(t,N,v)[[1]])
  )
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

#' @export
k_CAEX <- function(t=2,N,v=3){
  if (t>2 || !v==3) return(NA)
  k_TJcat(t,N,v)
}

#' @export
k_recBoseCA <- function(t=2, N, v, type="PCA"){
  if (t>2 || !v%in% primedat$q) return(NA)
  ## PCA: N = d*q^2 - (d-1)*q,   k = q^d + d*q^(d-1)
  ## CA:  N = d*q^2 - (d-1)*2,   k = (q+1)^d
  ##
  if (N < v^2) return(NA)
  if (type=="PCA"){
    ## get maximum d from N and q
    for (j in 1:10){
      if (N >= j*v^2 - (j-1)*v) d <- j else break
    }
    return(N_k_recursiveBose(v, d, type="PCA")["k"])
  }else ## type "CA"
  {
    ## get maximum d from N and q
    for (j in 1:10){
      if (N >= j*v^2 - (j-1)*2) d <- j else {
        break
      }
    }
    if (d>2) message("The assessment for type 'CA' is cautious\n, in some cases a smaller array might suffice")
    suppressMessages(return(N_k_recursiveBose(v, d, type="CA")["k"]))
  }
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

ks_CK_doubling <- function(t=3,N,v){
  if (!v %in% c(2,3)) return(NA)
  if(!t==3) return(NA)
  direct <- ks(3,N,v)
  khilf <- ceiling(max(ks)/2)
  if (v==3) minus <- suppressMessages(N_CAEX(k=khilf)) else
    minus <- N_KSK(khilf)
#### das ist noch gar nicht fertig !!!
  kdoubled <- ks(3, N-minus, v)
  kdoubled <- pmin(kdoubled, khilf)
  hilf <- Ns(3, khilf, v)
  Ndoubled <- hilf + minus
  doubled <- as.data.frame(rbind(hilf + plus))
  direct <- as.data.frame(rbind(Ns(t, k, v)))
  dplyr::bind_rows(direct=direct, doubled=doubled, .id="way")
}

## obtain sizes for product construction
## with large D1 and small D2
#' @export
ks_productCA <- function(t=2, N, v){
  ## the function compares direct ks for N
  ## with ks from subtracting the smallest "to-be-added" run number M-2
  ##   from N and multiplying the numbers of columns from the remaining
  ##   N - M + 2 with the number of columns l from the smallest CA for v
  if (!t==2) return(NA)
  M <- v^2
  if (v %in% primedat$q)
    l <- v+1
  else
    l <- 3 ## latin square
  Nhilf <- N - M + 2
  hilf <- ks(2, Nhilf, v)
  productCA <- as.data.frame(rbind(hilf*l))
  suppressMessages(direct <- as.data.frame(rbind(ks(t, N, v))))
  dplyr::bind_rows(direct=direct, productCA=productCA, .id="way")
}

## obtain sizes for product construction
## with large D1 and small D2
#' @export
ks_productCA <- function(t=2, N, v){
  ## the function compares direct ks for N
  ## with ks from subtracting the smallest "to-be-added" run number M-2
  ##   from N and multiplying the numbers of columns from the remaining
  ##   N - M + 2 with the number of columns l from the smallest CA for v
  if (!t==2) return(NA)
  M <- v^2
  if (v %in% primedat$q)
    l <- v+1
  else
    l <- 3 ## latin square
  Nhilf <- N - M + 2
  hilf <- ks(2, Nhilf, v)
  productCA <- as.data.frame(rbind(hilf*l))
  suppressMessages(direct <- as.data.frame(rbind(ks(t, N, v))))
  dplyr::bind_rows(direct=direct, productCA=productCA, .id="way")
}

