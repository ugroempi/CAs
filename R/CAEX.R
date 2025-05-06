#' CAEX: current optimal CAs of strength 2 with v=3
#'
#' Function for CAEX CAs of strength 2 with v=3.
#' These yield current optimal CAs.
#' The functions use stored CAs from Torres-Jimenez, Acevedo-Juarez
#' and Avila-George (2021), where necessary,
#' and product constructions from that paper, where possible.
#'
#' @rdname CAEX
#'
#' @aliases CAEX
#' @aliases k_CAEX
#'
#' @usage CAEX(k=NULL, N=NULL, t=2, v=3, maxk1=FALSE, ...)
#' @usage k_CAEX(N=NULL, t=2, v=3, ...)
#'
#' @param k positive integer or \code{NULL}; the requested number of factors;\cr
#' if both \code{k} and \code{N} are specified, \code{k} takes precedence;\cr
#' if both \code{k} and \code{N} are \code{NULL}, function \code{CAEX} throws an error.
#' @param N positive integer or \code{NULL}; the requested number of runs;
#' run sizes up to 50 are implemented;
#' these accommodate up to 112770(!) columns\cr
#' for \code{k_CAEX}, \code{N=NULL} implies that k for all possible
#' \code{N} (9, 11 to 50) \code{k} and \code{N} are returned (as part of
#' data.frame that also has \code{N} and the construction type).
#' @param t integer, the requested strength; currently,
#' only 2 is implemented
#' @param v integer, the requested number of levels; currently,
#' only 3 is implemented
#' @param maxk1 logical: should it be attempted to achieve maximum k1 for
#' the output's PCA structure? \code{TRUE} may substantially increase run time.
#' @param ... currently not used
#'
#' @section Details:
#' Function \code{CAEX} provides only strength 2 CAS, and currently only for 3-level columns.
#' These are the current best ones, i.e., the necessary run size for a desired number of
#' columns \code{k} can be obtained from function \code{\link{eCAN}},
#' or from function \code{\link{N_TJcat}}.
#'
#' Function \code{CAEX} implements the arrays from Torres-Jimenez, Acevedo-Juarez and
#' Avila-George (2021). It uses internally-stored CAs (in object \code{CAEX_CAs}),
#' where the last construction step is based on numerical column extension
#' (N=11 to 33 and 35), and product constructions as indicated in
#' Table 4 of Torres-Jimenez, Acevedo-Juarez and Avila-George (2021)
#' for N = 34 and 36 to 50. The information on how to proceed comes from the internal
#' list \code{CAEX_lineages}. In all cases, the returned strength 2 (P)CAs have the same
#' dimensions as those listed in the paper; where products were calculated by function
#' \code{CAEX}, the \code{k1} is per default smaller than that provided in the paper,
#' as it is not undertaken to optimize it beyond moving the already suitable columns
#' to the left.
#'
#' If the designs are intended for use in further PCAxPCA constructions, maximization
#' of the k1 value in the resulting PCA structure is desirable.
#' Setting \code{maxk1} to TRUE will attempt to achieve
#' the k1 indicated in Torres-Jimenez et al. (2021).
#' This is attempted by using the internal function \code{find_clique}, which
#' searches for a subset of v rows with at least k1 columns with
#' distinct levels. If successful, these rows are moved to the top,
#' and function \code{CA_to_PCA} yields the split k = k1+k2 with k1 as stated in the paper.\cr
#' The current implementation is not ideal regarding flexible values (coded as \code{NA})
#' among the possibilities for increasing k1 (a single NA in a v-tuple can be made use of,
#' two or more NAs are ignored, due to R's handling of missing values in assessing uniqueness);
#' There are not many flexible values in the stored ingredient CAs
#' (none for N=11to 15, 17, 18, 20 and 21, 1 for N=23, 24 or 27, 2 for N=16, 4 for N=26,
#' 11 for N=19 or 29, 12 for N=22, 17 for N=30, 48 for N=25 or 32, 50 for N=33, 65 for N=31,
#' 69 for N=35, 197 for N=28). It is thus conceivable, if not likely, that it may be impossible to achieve the
#' k1 value attained in Torres-Jimenez et al. (2021); in that case, the best possible attainable
#' k1 value will be produced, which is likely larger than the default, but may require a long run time.
#' The search may take a long time, if a successful v-tuple of rows is found late in the
#' search path over the \code{choose(N, v)} tuples. Thus, this search should only be triggered
#' if a large \code{k1} is really needed; it is not relevant for the coverage properties of
#' the design itself.
#'
#' @returns Function \code{CAEX} returns a CA (matrix of class \code{ca})
#' with levels coded from 0 to v-1, and flexible positions (if any) denoted as \code{NA}.
#' The returned object has an attribute \code{PCAstatus} that shows the PCA status,
#' with k1 and k2 obtained from applying function \code{\link{CA_to_PCA}}
#' to the result in case of product calculation.
#' The splits k=k1+k2 are worse than the splits reported in Torres-Jimenez et al. (2021)
#' for the product constructions (see also listing in the Examples section),
#' unless \code{maxk1=TRUE} is chosen, at the expense of much longer run time.\cr
#' Function \code{k_CAEX} returns a data frame with columns \code{N},
#' \code{k} and \code{method},
#' where \code{k} holds the maximum number of columns for run sizes up to \code{N},
#' and \code{method} indicates, whether the design is the \code{Bose} array (\code{N=9}),
#' a stored array, or is obtained from stored arrays via a \code{CAxCA} or a \code{PCAxPCA}
#' construction. For \code{N=NULL}, a \code{data.frame}
#' of all CAEX constructions is returned.
#'
#' @references Torres-Jimenez, Acevedo-Juarez and Avila-George (2021)
#'
#' @examples
#' ## a stored CA
#' D <- CAEX(12) ## k=12 columns
#' eCAN(2,12,3)
#' coverage(D, 2)
#'
#' ############################################################
#' ## CAs obtained using function productCA as the workhorse
#'
#' ## from two CAs without flexible values
#' D <- CAEX(N=34)
#' attributes(D)
#' is.PCA(D)
#'
#' ## from two CAs with flexible values
#' ## not run because of slightly longish run time
#' \dontrun{
#' D <- CAEX(N=43)
#' attributes(D)
#' }
#' ############################################################
#'
#' ##############################################################
#' ## available constructions from function CAEX
#' ##############################################################
#' \dontrun{
#' # this would run too long for some arrays
#' # longest about 07 min for N=49 (using productCA, like two other cases)
#' # about 13 s for N=50 (using productPCA, like most cases)
#' allCAEX <- lapply(11:50, function(obj) CAEX(N=obj))
#'
#' dimsAchieved <- t(sapply(allCAEX, dim))
#' colnames(dimsAchieved) <- c("N","k")
#' k1Achieved <- sapply(allCAEX, function(obj) attr(obj, "PCAstatus")$k1)
#' (CAEXinfos <- cbind(as.data.frame(cbind(dimsAchieved, k1=k1Achieved)),
#'      method=unlist(lapply(CAs:::CAEX_lineages[[1]][[1]],
#'            function(obj) ifelse(is.character(obj), "stored", obj[1]))),
#'      k1_paper=c(4,4,6,7,14,15,22,31,34,45,69,76,99,144,138,185,294,295,414,
#'                 582,655,872,1327,1426,1938,2918,3210,4163,6262,6842,8838,
#'                 12216,13485,18894,27606,30271,39387,55151,60884,76842)
#'                 ))
#' ## see optimality by comparing with
#' ## TJcat[which(TJcat[,"t"]==2 & TJcat[,"v"]==3 & TJcat[,"N"]<=50)[-1],c("k","N")])
#'
#'      N      k    k1  method k1_paper
#' # 11 11      5     4  stored        4
#' # 12 12      7     4  stored        4
#' # 13 13      9     6  stored        6
#' # 14 14     10     7  stored        7
#' # 15 15     20    14  stored       14
#' # 16 16     21    15  stored       15
#' # 17 17     29    22  stored       22
#' # 18 18     46    31  stored       31
#' # 19 19     49    34  stored       34
#' # 20 20     63    45  stored       45
#' # 21 21     93    69  stored       69
#' # 22 22    107    76  stored       76
#' # 23 23    138    99  stored       99
#' # 24 24    199   144  stored      144
#' # 25 25    216   138  stored      138
#' # 26 26    288   185  stored      185
#' # 27 27    435   294  stored      294
#' # 28 28    449   295  stored      295
#' # 29 29    610   414  stored      414
#' # 30 30    878   582  stored      582
#' # 31 31    964   655  stored      655
#' # 32 32   1308   872  stored      872
#' # 33 33   1964  1327  stored     1327
#' # 34 34   2116   460   CAxCA     1426  ## < 1 second,
#'                                        ## or 95 seconds with maxk1=TRUE
#' # 35 35   2700  1938  stored     1938
#' # 36 36   3918  2139 PCAxPCA     2918  ## < 1 second,
#'                                        ## or > 8 minutes with maxk1=TRUE
#' # 37 37   4457  2356 PCAxPCA     3210
#' # 38 38   5763  3069 PCAxPCA     4163
#' # 39 39   8329  4464 PCAxPCA     6262
#' # 40 40   9207  5244 PCAxPCA     6842
#' # 41 41  11898  6831 PCAxPCA     8838
#' # 42 42  17895  9114 PCAxPCA    12216  ## approx 1.3 seconds (productPCA),
#'                                        ##   or approx. 10.5 min with maxk1=TRUE
#' # 43 43  20010  4350   CAxCA    13485  ## approx 16 seconds (productCA slower),
#'                                        ##   or approx 11.5 min with maxk1=TRUE
#' # 44 44  25317 14256 PCAxPCA    18894  ## approx 1.9 seconds (productPCA)
#'                                        ##   or more than 40 min with maxk1=TRUE
#' # 45 45  37071 20286 PCAxPCA    27606
#' # 46 46  42174 22344 PCAxPCA    30271
#' # 47 47  54531 29106 PCAxPCA    39387
#' # 48 48  80789 41137 PCAxPCA    55151
#' # 49 49  90344 19640   CAxCA    60884  ## approx 7 min (productCA slower)
#' # 50 50 112770 60078 PCAxPCA    76842  ## approx 13 seconds (productPCA)
#'
#' }

#' @export
CAEX <- function(k=NULL, N=NULL, t=2, v=3, maxk1=FALSE, ...){
  call <- sys.call()
  stopifnot(t==2, v==3)
  if (is.null(k) && is.null(N)) stop("at least one of k and N must be specified")
  stopifnot(is.logical(maxk1))

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
    if (N==10){
      N <- 9
      message("N has been reduced to 9, the largest possible N <= 10")
    }
    if (N>50){
      N <- 50
      message("N has been reduced to 50, the current maximum for CAEX.")
    }
    k <- min(TJcat[which(TJcat[,"t"]==2 & TJcat[,"v"]==3 & TJcat[,"N"]>=N), "k"])
  }
  if (t==2 && v==3){
    ## current only case, precaution for future extensions
    ## arrays only
    if (N==v^2 || k==v+1) D <- SCA_Bose(v) else{
      lineages <- CAEX_lineages[[as.character(t)]][[as.character(v)]]
      lineage <- lineages[[as.character(N)]]
      if (is.character(lineage)){
        D <- CAEX_CAs[[lineage]]
      }else{
        if (maxk1) k1_paper <- lineage$k1_paper
        D1 <- CAEX_CAs[[lineage$one]]
        ## take care of the case where D1 is in turn constructed
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
            D <- productCA(D1, D2, one_is_enough=TRUE)
        ## for t=2 and v=3, one_is_enough is speedier and does not deteriorate N
        ## make k1 as large as possible in the k1 to k2 split
        ## by swapping design columns

        ## prior step of picking rows with suitable clique size,
        ## which is known from the paper
        if (maxk1){
          hilf <- find_clique(D, target=k1_paper)
          ## treat the case for which NA values occur
          ## in the first few rows of the D to hand to CA_to_PCA
          if (any(is.na(D[hilf,]))){
            NAcols <- which(apply(D[hilf,], 2, function(obj) any(is.na(obj))))
            for (c in NAcols){
              ## treat columns that will be moved to the front and have levels
              ## swapped by CA_to_PCA
              if (length(unique(D[hilf,c]))==v){
                ## at most one NA, because all NA are treated as identical
                NApos <- which(is.na(D[hilf,c]))    ## only one value possible
                D[NApos,c] <- setdiff(0:(v-1), D[hilf,c])  ## the single missing level
              }
            }
          }
          D <- D[c(hilf, setdiff(1:N, hilf)),]
        }
        D <- CA_to_PCA(D)
      }
    }
  } ## end of t=2, v=3
  if (k < ncol(D)){
    korig <- ncol(D)
    D <- D[,1:k]
    if (!is.null(attr(D, "PCAstatus"))){
      hilf <- attr(D, "PCAstatus")
      if (hilf$type %in% c("PCA", "SCA")){
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

#' @export
k_CAEX <- function(N=NULL, t=2, v=3, ...){
  ## prepare
  states <- cbind(k=c(4,5,7,9,10,20,21,29,46,49,63,93,107,138,199,216,288,435,449,610,878,964,1308,1964,2116,2700,
                      3918,4457,5763,8329,9207,11898,17895,20010,25317,37071,42174,54531,80789,90344,112770),
                  N=c(9,11:50))
  constrs <- c("Bose", rep("stored", 23), "CAxCA", "stored", rep("PCAxPCA", 7),
               "CAxCA", rep("PCAxPCA", 5),
               "CAxCA", "PCAxPCA")
  df <- data.frame(states[,2:1], method=constrs)
  ## all
  if (is.null(N)) return(df) else{
    ## one only
    if (length(N)>1){
      message("used first element of N")
      message("use N=NULL for seeing all construction sizes")
      N <- N[1]
    }
    stopifnot(is.numeric(N))
    stopifnot(N>=9)
    if (N==10){
      message("N reduced to 9 (largest N <=10 for the construction")
      N <- 9
    }
    if (N>50){
      message("N reduced to 50 (largest implemented N for the construction")
      N <- 50
    }
    N <- floor(N)
    constr <- CAEX_lineages[[as.character(t)]][[as.character(v)]][[as.character(N)]]
    if (is.character(constr)) {
      kchar <- as.numeric(rev(strsplit(constr,".",fixed=TRUE)[[1]])[1])
      k <- states[,"k"][which(states[,"N"]==N)]
      stopifnot(kchar==k) ## programming mistake
      return(df[which(df$k==k),])
    }else{
      k <- states[,"k"][which(states[,"N"]==N)]
      storedmethod <- df[which(df[,"N"]==N),"method"]
      method <- constr$method
      stopifnot(method==storedmethod)
      chr1 <- constr$one
      chr2 <- constr$two
      return(data.frame(N=N, k=k, method=method, one=chr1, two=chr2))
    }
  }
}

find_clique <- function(D, target, v=NULL, ...){
  ## it is assumed that D is uniform with levels 0 to v-1
  ## with potential flexible values coded as NA
  N <- nrow(D); k <- ncol(D)
  if (is.null(v)) v <- max(D, na.rm=TRUE) + 1
  stopifnot(v > 1)
  sets <- nchoosek(N, v)
  ## would this be better?
  ## so far did not seem so.
  fun <- function(x){
    length(unique(setdiff(x, NA))) + sum(is.na(x)) == v
  }
  best <- 1:v
  k1best <- sum(apply(D[1:v,], 2, fun)) #function(obj) length(unique(obj))==v))
  for (i in 2:ncol(sets)){
    cur <- sum(apply(D[sets[,i],], 2, fun)) #function(obj) length(unique(obj))==v))
    if (cur > k1best) {
      best <- sets[,i]
      k1best <- cur
    }
    if (k1best >= target) break
  }
  best
}
