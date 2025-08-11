### Construction of Lobb et al. 2012 (LCDST) based on Z_v-1
### and circular matrix created from a cover starter vector
###    taken fom LCDSTStarters
#'
#' Strength 2 CAs constructed by cover starters (LCDST)
#'
#' Function to construct a strength 2 CA based on cover starters according to Lobb et al. 2012.
#' k and v must be specified, the number of runs is ((v-f)*k+CAN_(2,k,f),
#' where f is the number of fixed points (assuming that the CAN array is available,
#' which is the case for f=1,2,3 (eCAN in case of f=3))
#'
#' @rdname CS_LCDST
#'
#' @aliases CS_LCDST
#' @aliases N_CS_LCDST
#' @aliases k_CS_LCDST
#'
#' @usage CS_LCDST(k, v, start0=TRUE, starter=NULL, ...)
#' @usage N_CS_LCDST(t=2, k, v)
#' @usage k_CS_LCDST(t=2, N, v)
#'
#' @param k number of factors
#' @param v number of levels, \code{v} >= 3
#' @param start0 should the output design's levels start with zero (default) or 1?
#' @param starter an optional cover starter vector that starts with a fixed value, denoted by a number larger than 100.\cr
#' This may be useful for cases with k larger than implemented in the stored cover starters.\cr
#' CAUTION: Responsibility for achieving the coverage rests with the user, please check!
#' @param t requested strength (must be 2, added of unified interface with other N_ and k_ functions)
#' @param N affordable number of runs
#' @param ... currently not used
#'
#' @returns \code{CS_LCDST} returns the smallest possible strength 2 CA from the Lobb et al. (2012) construction for \code{k} experimental factors in \code{v} levels; it is an \code{N x k} matrix with entries from 0 to \code{v-1} or from 1 to \code{v}, depending on \code{start0}.\cr
#' \code{N_CS_LCDST} returns the minimum run size needed for the LCDST construction with \code{k} factors at \code{v} levels, or \code{NA}, if no such construction exists.\cr
#' \code{k_CS_LCDST} returns the maximum number of \code{v}-level factors that can be accommodated in at most \code{N} runs, or \code{NA}, if no such construction exists.
#'
#' @section Details:
#' Lobb et al. (2012) denote designs with the factors as the columns and the runs as the rows,
#' like in this package, and contrary to Meagher and Stevens (2005).
#'
#' The construction uses starter vectors of \code{k} elements with first element a fixed position,
#' denoted as a number >100, depending on \code{v}: the fixed points are taken as the last f levels,
#' and they will be made valid levels \code{v-1} to \code{v-f} eventually by subtracting 100.
#'
#' The starter vector is cycled into a \code{k x k} matrix \code{M} (internal function \code{circular}),
#' which is one key ingredient of the construction. The other key ingredient is the additive group \code{Z} of
#' \code{v-f} elements, which is repeatedly applied to the non-fixed elements.
#' The resulting matrix consists in vertical concatenations of those \code{v-f} matrices,
#' augmented by the smallest possible CA for \code{k} f-level columns in order to cover
#' interactions among the fixed levels (a constant row with the only fixed level in case of f=1).
#'
#' Functions \code{N_CS_LCDST} and \code{k_CS_LCDST} return the minimum \code{N} or the maximum
#' \code{k} for the requested scenario. In case of \code{N_CS_LCDST}, the default to use the
#' available starters can be overridden by the argument \code{theoretical}, which allows
#' to inspect how many runs can be obtained by a suitable user-specified starter vector
#' with the perfect additional ingredient.
#'
#' @seealso [LCDSTCombis]
#' @seealso [LCDSTStarters]
#'
#' @references Lobb et al. (2012)
#'
#' @examples
#' ## example for k=5, v=3, f=1
#' starter <- c(102, 0, 1, 1, 0)
#' coverage(CS_LCDST(5, 3, starter=starter), 2)
#'
#' ## example for k=20, v=8, f=3
#' ## levels 0 to 7, with levels 5 to 7 fixed
#' ## array used by function CS_LCSDT for interactions among fixed points
#' addfix <- CAEX(20)  ##  20 3-level columns at strength 2
#' coverage(addfix,2)
#' dim(addfix)         ## 15 rows
#' CA115.2.20.8 <- CS_LCDST(20, 8) ## 115 runs, 20 columns at 8 levels
#' dim(CA115.2.20.8)   ## (v-f)*k + nrow(addfix) = (8 - 3)*20 + 15 rows
#' coverage(CA115.2.20.8,2)
#'
#' N_CS_LCDST(k=20, v=8)  ## 115 runs are needed for 20 8-level columns
#' k_CS_LCDST(N=120, v=8) ## with up to 120 runs, 20 8-level columns are possible
#' eCAN(2, 20, 8) ## post-processing with simulated annealing may loose seven rows
#'
#' ## f larger than 3
#' ## v=11 levels, f=4, k=26
#' D <- CS_LCDST(26, 11)
#' dim(D)
#' N_CS_LCDST(2,26,11)
#' LCDSTCombis[which(LCDSTCombis$k==26 & LCDSTCombis$v==11),]
#'

#' @export
CS_LCDST <- function(k, v, start0=TRUE, starter=NULL, ...){
  Call <- sys.call()
  ## checks
  stopifnot(is.numeric(k), is.numeric(v))
  if (v == 2) stop(paste0("For v=2, use function KSK .\n It yields a strength 2 CA in ",
                          N_KSK(k), "runs."))
  stopifnot(v >= 3)
  stopifnot(v %% 1 == 0)
  stopifnot(k > 2)
  stopifnot(is.logical(start0))
  ## a specified starter
  if (!is.null(starter)){
    stopifnot(length(starter) >= k)
    stopifnot(starter[1] > 100)
  }
  if (is.null(starter)){
    N <- N_CS_LCDST(t=2, k, v)
    if (is.na(N)) stop("This combination of k and v cannot be accommodated.")
    k_lookup <- k_CS_LCDST(t=2,N,v)
    if (is.na(k_lookup)) stop("Unexpected error.")
    ## retrieve starter
    vc <- as.character(v)
    kc <- as.character(k_lookup)

    starter <- LCDSTStarters[[vc]][[kc]]
    if (is.null(starter)) stop("not implemented")
  }
  nNA <- sum(is.na(starter))
  f <- length(unique(starter[which(starter>100)])) ## fixed values
  ## there is no case with more NA values than unfixed levels
  if (nNA>0) starter[is.na(starter)] <- (0:(v-f-1))[1:nNA]

  ## create array
  G <- createG_LCDST(v, f)
  aus <- createCS_LCDST(starter, v, G)[,1:k]
  if (!start0) aus <- aus + 1
  attr(aus, "origin") <- "CS_LCDST"
  attr(aus, "t") <- 2
  attr(aus, "f") <- f
  attr(aus, "Call") <- Call
  if (f>1)
  hilf <- Ns(2,k,f)
  if (f>3){
    attr(aus, "method_D-for-f") <- names(hilf)[
                             which.min(hilf[-length(hilf)])]
    attr(aus, "Date") <- Sys.Date()
    attr(aus, "CAs version") <- packageVersion("CAs")
  }
  aus
}

## the N and k functions
#' @export
N_CS_LCDST <- function(t=2, k, v){
  if (!t==2 || v==2) return(NA)
  if (v == 2) stop(paste0("For v=2, use function KSK .\n It yields the optimum run size, which is ", N_KSK(k), "."))
  ## checks
  stopifnot(is.numeric(k), is.numeric(v))
  stopifnot(k%%1==0, v%%1==0)
  stopifnot(k>=2, v>=3)
  hilf <- LCDSTCombis[LCDSTCombis[,"v"]==v &
                      LCDSTCombis[,"k"]>=k,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  # would be a solution if N were eCAN-based
  # if (hilf[1,"f"]==4) return(hilf[1,"N"] + bestN(2,k,4)-eCAN(2,k,4)$CAN)
  hilf[1,"N"]
}

#' @export
k_CS_LCDST <- function(t=2, N, v){
  if (!t==2) return(NA)
  stopifnot(is.numeric(N), is.numeric(v))
  if (v == 2) stop(paste0("For v=2, use function KSK .\n It yields the maximum number of factors, which is ", k_KSK(N), "."))
  hilf <- LCDSTCombis[LCDSTCombis[,"v"]==v & LCDSTCombis[,"N"]<=N,,drop=FALSE]
  # would be a solution if N were eCAN-based
  # if (any(hilf[,"f"]==4)){
  #   hilf[which(hilf[,"f"]==4),"N"] <- hilf[which(hilf[,"f"]==4),"N"] +
  #     sapply(hilf[which(hilf[,"f"]==4),"k"],
  #            function(obj) bestN(2, obj, 4) - eCAN(2, obj, 4)$CAN)
  #   hilf <- hilf[which(hilf[,"N"] <=N),]
  # }
  if (nrow(hilf)==0) return(NA)
  posmax <- which.max(hilf[,"k"])
  hilf[posmax,"k"]
}

### auxiliary functions

createG_LCDST <- function(v, f){
  ## creates function that uses the additive group
  ## except for the fixed elements
  ## levels are coded from 0 to v-f-1
  function(x, f) ifelse(x>100, x, (x+1) %% (v-f))
}

createCS_LCDST <- function(starter, v, G, ...){
  # worker function, checks happen in the interface function
  # returns object in 1 to v coding
  # starter: starter vector of length k with values from 0 to v-1, first value 0,
  #          retrieved from stored list of lists MeagherStevensStarters
  # v: number of levels
  # G: a permutation function, created with createG
  # ...: currently not used
  #
  ## create M by cycling the starter vector
  ## 1 to v coding of levels
  M <- circular(starter)
  k <- length(starter)
  f <- length(unique(starter[which(starter>100)]))

  ## then M
  ## then M with permuted non-fixed elements (by G)
  mat <- M
  for (i in 1:(v-f-1)){
    M <- G(M, f)
    mat <- rbind(mat, M)
  }
  ## now interactions among fixed points
  if (f==1) mat <- rbind(mat, v-1)
  if (f==2) mat <- rbind(mat, KSK(k=k) + v - f)
  if (f==3) mat <- rbind(mat, CAEX(k=k) + v - f)
  if (f>3) mat <- rbind(mat, bestCA(2,k,f) + v - f)

  fixed <- which(mat>100)
  mat[fixed] <- mat[fixed]-100  ## fix fixed values
  mat
}

## internal function for searching starters for small cases
## has to be edited for adapting it to the search in question
## was only ever used manually
search_CS <- function(v){
  ## k is implemented through the length of the starter
  for (i1 in 1:(v-1))
    for (i2 in 1:(v-1))
      for (i3 in 1:(v-1))
        for (i4 in 1:(v-1))
          for (i5 in 1:(v-1)){
    starter <- c(0,1,1,1,1,i1,i2,i3,i4,i5)
    k <- length(starter)
    M <- circular(starter) + 1
    mat <- cbind(1,M)
    G <- createG(v)
    for (i in 2:(v-1)){
      M <- matrix(G(M), nrow=k)
      mat <- cbind(mat, M)
    }
    ## transpose, because rows and columns swapped versus paper
    ## the dimension is correct, but the coverage is not
    hilf <- coverage(t(mat)-1, 2)
    if (all(unlist(hilf)==1)) {
      print(starter)
      return(list(v=v, k=k, starter=starter, G=G))
    }
          }
  message("The search for k =", k, "columns in v =", v, "levels was not successful.")
}
