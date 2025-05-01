#' CAEX: current optimal CAs of strength 2 with v=3
#'
#' Functions for CAEX CAs of strength 2 with v=3.
#' These yield current optimal CAs.
#' The functions use stored CAs from Torres-Jimenez, Acevedo-Juarez
#' and Avila-George (2021), where necessary,
#' and product constructions from that paper, where possible.
#'
#' @rdname CAEX
#'
#' @aliases CAEX
#'
#' @usage CAEX(k=NULL, N=NULL, t=2, v=3, ...)
#'
#' @param k positive integer or \code{NULL}; the requested number of factors;\cr
#' if both \code{k} and \code{N} are specified, \code{k} takes precedence;\cr
#' if both \code{k} and \code{N} are \code{NULL}, function \code{CAEX} throws an error.
#' @param N positive integer or \code{NULL}; the requested number of runs; run sizes up to 50 are implemented;
#' these accommodate up to 112770(!) columns\cr
#' @param t integer, the requested strength; currently, only 2 is implemented
#' @param v integer, the requested number of levels; currently, only 3 is implemented
#' @param ... currently not used
#'
#' @section Details:
#' Function \code{CAEX} implements the arrays from Torres-Jimenez, Acevedo-Juarez and Avila-George (2021).
#' It uses internally-stored CAs (in object \code{CAEX_CAs}) where the last construction step is based on numerical column extension
#' (N=11 to 33 and 35), and product constructions as indicated in Table&nbsp; of Torres-Jimenez,
#' Acevedo-Juarez and Avila-George (2021) for N = 34 and 36 to 50.
#' The information on how to proceed comes from the internal list \code{CAEX_lineages}.\cr
#' The arrays have values starting with 0; flexible values, if any, are denoted as \code{NA}.
#' Function \code{CAEX} provides only strength 2 CAS, and currently only for 3-level columns.
#' These are the current best ones, i.e., the nessecary run size for a desired number of columns \code{k} can be obtained from function \code{\link{eCAN}},
#' or from function \code{\link{N_TJcat}}.
#'
#' @returns Function \code{CAEX} returns a CA with levels coded from 0 to v-1, and flexible positions
#' denoted as \code{NA} and an attribute that shows the PCA status (with k1 and k2 arising from the product
#' without postprocessing in case of product constructions).\cr
#'
#' @references Torres-Jimenez, Acevedo-Juarez and Avila-George (2021)
#'
#' @examples
#' ## a stored CA
#' CAEX(12) ## k=12 columns
#' eCAN(2,12,3)
#'
#' ## a CA obtained from other CAs without flexible values by function \code{\link{productCA1}}
#' D <- CAEX(N=34) ## should have N=34 runs but currently has 35
#'                 ## --> needs improvement of function productCA1
#' attributes(D)
#' is.PCA(D)
#'
#' ## a CA obtained from other CAs with flexible values by function \code{\link{productCA1}}
#' ## not run because of slightly longish run time
#' \dontrun{
#' D <- CAEX(N=43) ## should have N=43 runs, but has 44
#'                 ## improvement to productCA1 needed
#' attributes(D)
#' }
#'
# # this would run a long time for the larger arrays
# # about 10 min for the largest one
# dimsAchieved <- sapply(11:50, function(obj) dim(CAEX(N=obj)))
# cbind(dimsAchieved, TJcat[which(TJcat[,"t"]==2)[-1],]))
#                             t      k v  N
# ca.2.3^5.txt      11      5 2      5 3 11
# ca.2.3^7.txt      12      7 2      7 3 12
# ca.2.3^9.txt      13      9 2      9 3 13
# ca.2.3^10.txt     14     10 2     10 3 14
# ca.2.3^20.txt     15     20 2     20 3 15
# ca.2.3^21.txt     16     21 2     21 3 16
# ca.2.3^29.txt     17     29 2     29 3 17
# ca.2.3^46.txt     18     46 2     46 3 18
# ca.2.3^49.txt     19     49 2     49 3 19
# ca.2.3^63.txt     20     63 2     63 3 20
# ca.2.3^93.txt     21     93 2     93 3 21
# ca.2.3^107.txt    22    107 2    107 3 22
# ca.2.3^138.txt    23    138 2    138 3 23
# ca.2.3^199.txt    24    199 2    199 3 24
# ca.2.3^216.txt    25    216 2    216 3 25
# ca.2.3^288.txt    26    288 2    288 3 26
# ca.2.3^435.txt    27    435 2    435 3 27
# ca.2.3^449.txt    28    449 2    449 3 28
# ca.2.3^610.txt    29    610 2    610 3 29
# ca.2.3^878.txt    30    878 2    878 3 30
# ca.2.3^964.txt    31    964 2    964 3 31
# ca.2.3^1308.txt   32   1308 2   1308 3 32
# ca.2.3^1964.txt   33   1964 2   1964 3 33
# ca.2.3^2116.txt   35   2116 2   2116 3 34  *** productCA1
# ca.2.3^2700.txt   35   2700 2   2700 3 35
# ca.2.3^3918.txt   36   3918 2   3918 3 36
# ca.2.3^4457.txt   37   4457 2   4457 3 37
# ca.2.3^5763.txt   38   5763 2   5763 3 38
# ca.2.3^8329.txt   39   8329 2   8329 3 39
# ca.2.3^9207.txt   40   9207 2   9207 3 40
# ca.2.3^11898.txt  41  11898 2  11898 3 41
# ca.2.3^17895.txt  42  17895 2  17895 3 42
# ca.2.3^20010.txt  44  20010 2  20010 3 43  *** productCA1
# ca.2.3^25317.txt  44  25317 2  25317 3 44
# ca.2.3^37071.txt  45  37071 2  37071 3 45
# ca.2.3^42174.txt  46  42174 2  42174 3 46
# ca.2.3^54531.txt  47  54531 2  54531 3 47
# ca.2.3^80789.txt  48  80789 2  80789 3 48
# ca.2.3^90344.txt  50  90344 2  90344 3 49  *** productCA1
# ca.2.3^112770.txt 50 112770 2 112770 3 50
#
# rownames(dimsAchieved) <- c("N_achieved", "k_achieved")
# dimerg <- cbind(t(dimsAchieved), TJcat[which(TJcat[,"t"]==2)[-1],])
#

#' @export
CAEX <- function(k=NULL, N=NULL, t=2, v=3, ...){
  call <- sys.call()
  stopifnot(t==2, v==3)
  if (is.null(k) && is.null(N)) stop("at least one of k and N must be specified")
  if (!is.null(k) && !is.null(N)){
    N <- NULL
    message("If both k and N are specified, k takes precedence,\n and the minimum N that accommodates k factors is chosen.")
  }
  if (is.null(N)){
    ## check input k
      stopifnot(length(k)==1)
      stopifnot(is.numeric(k))
      stopifnot(k%%1==0)
      stopifnot(k>=1)
    ## determine N
      N <- N_TJcat(t, k, v)
      if (is.na(N)) stop("The requested k =", k, " cannot be accommodated.")
  }else{
    ## check input N
    stopifnot(length(N)==1)
    stopifnot(is.numeric(N))
    stopifnot(N%%1==0)
    stopifnot(N>=9)
    if (N>50){
      N <- 50
      message("N has been reduced to 50, the current maximum for CAEX.")
    }
    k <- min(TJcat[which(TJcat[,"t"]==2 & TJcat[,"v"]==3 & TJcat[,"N"]>=N), "k"])
  }
  if (t==2 && v==3){
    ## current only case, precaution for future extensions
    ## arrays only
    if (N==v^2) return(SCA_Bose(v))
    lineages <- CAEX_lineages[[as.character(t)]][[as.character(v)]]
    lineage <- lineages[[as.character(N)]]
    if (is.character(lineage)){
      D <- CAEX_CAs[[lineage]]
    }else{
      D1 <- CAEX_CAs[[lineage$one]]
      ## take care of the case where D1 is in turn constructe
      ## at present not relevant (30 April 2025)
      if (is.null(D1)) D1 <- CAEX(k=as.numeric(strsplit(lineage$one, ".",
                                                        fixed=TRUE)[[1]][[4]]))
      D2 <- CAEX_CAs[[lineage$two]]
      ## take care of the case where D2 is in turn constructe
      ## at present not relevant (30 April 2025)
      if (is.null(D2)) D2 <- CAEX(k=as.numeric(strsplit(lineage$two, ".",
                                                        fixed=TRUE)[[1]][[4]]))
      if (lineage$method=="PCAxPCA")
        D <- productPCA(D1, D2) else
          D <- productCA1(D1, D2)
      ## make k1 as large as possible in the k1 to k2 split
      ## by swapping design columns
      D <- CA_to_PCA(D)
    }
  } ## end of t=2, v=3
  if (k < ncol(D)){
    korig <- ncol(D)
     D <- D[,1:k]
     if (!is.null(attr(D, "PCAstatus"))){
       hilf <- attr(D, "PCAstatus")
       if (hilf$PCAstatus$type %in% c("PCA", "SCA")){
         hilf$k2 <- max(0, hilf$k2 - (korig - k))
         hilf$k1 <- min(hilf$k1, k)
       }
     attr(D, "PCAstatus") <- hilf
     }
  }
  if (k < ncol(D)) D <- D[,1:k]
  hilf <- is.PCA(D, flexible=NA)
    if (hilf) attr(D, "PCAstatus") <- attr(hilf, "PCAstatus")
  if (!is.null(attr(hilf, "flexible")))
    attr(D, "flexible") <- attr(hilf, "flexible")
  class(D) <- c("ca", class(D))
  attr(D, "Call") <- call
  attr(D, "lineage") <- lineage
  D
}

