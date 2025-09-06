### Construction of Colbourn et al. 2006 (CMMSSY) based on Z_v-1
### and circular matrix created from a cover starter vector
###    taken fom CMMSSYStarters
###
### This also has some starters from LCDST
### However, it is easier - therefore it could be good to change Meagher and Stevens to
### this approach with Inf, because it would avoid using package permutations
###
#' Strength 2 CAs constructed by cover starters (CMMSSY)
#'
#' Function to construct a strength 2 CA based on cover starters according to Colbourn et al. 2006.
#' k and v must be specified, the number of runs is ((v-1)*k+1) (possibly for a larger value of k
#' if the construction does not exist for the requested k); recursive constructions based on such
#' arrays are also included
#'
#' @rdname CS_CMMSSY
#'
#' @aliases CS_CMMSSY
#' @aliases N_CS_CMMSSY
#' @aliases k_CS_CMMSSY
#'
#' @usage CS_CMMSSY(k, v, start0=TRUE, starter=NULL, ...)
#' @usage N_CS_CMMSSY(t=2, k, v, theoretical=FALSE)
#' @usage k_CS_CMMSSY(N, v)
#'
#' @param k number of factors
#' @param v number of levels, \code{v} >= 3
#' @param start0 should the output design's levels start with zero (default) or 1?
#' @param starter an optional cover starter vector that starts with \code{Inf}.\cr
#' This may be useful for cases with k larger than implemented in the stored cover starters.\cr
#' CAUTION: Responsibility for achieving the coverage rests with the user, please check!
#' @param t requested strength (has to be 2)
#' @param theoretical logical: default requests size available from using stored cover starters,\cr
#' change to \code{TRUE} for the size assuming existence of a cover starter for \code{k} and \code{v}
#' @param N affordable number of runs
#' @param ... currently not used
#'
#' @returns \code{CS_MMSSY} returns the smallest possible strength 2 CA from the Colbourn et al. (2006) construction for \code{k} experimental factors in \code{v} levels; it is an \code{N x k} matrix with entries from 0 to \code{v-1} or from 1 to \code{v}, depending on \code{start0}.\cr
#' \code{N_CS_MMSSY} returns the minimum run size needed for the CMMSSY construction with \code{k} factors at \code{v} levels, or \code{NA}, if no such construction exists.\cr
#' \code{k_CS_MMSSY} returns the maximum number of \code{v}-level factors that can be accommodated in at most \code{N} runs, or \code{NA}, if no such construction exists.
#'
#' @section Details:
#' Colbourn et al. (2006) denote designs with the factors as the columns and the runs as the rows, like in this package, and contrary to Meagher and Stevens (2005).
#'
#' The construction uses a distinct starter vectors of \code{k} elements with first element \code{Inf} (fixed position)
#' and subsequent elements from 1 to \code{v-1}.
#'
#' The starter vector is cycled into a \code{k x k} matrix \code{M} (internal function \code{circular}),
#' which is one key ingredient of the construction. The other key ingredient is addition modulo \code{v-1}.
#' Finally, the infinite element is set to \code{v-1}.\cr
#'
#' Functions \code{N_CS_CMMSSY} and \code{k_CS_MMSSY} return the minimum \code{N} or the maximum
#' \code{k} for the requested scenario. In case of \code{N_CS_MMSSY}, the default to use the
#' available starters can be overridden by the argument \code{theoretical}, which allows
#' to inspect how many runs can be obtained by a suitable user-specified starter vector.
#'
#' Besides the starters from Colbourn et al. (2006), the function also returns arrays that
#' are obtained by combining pairs, triples or quadruples of such arrays,
#' using function \code{\link{productPCA}}, once (for pairs) or repeatedly.
#' It would be possible to drive this further, but the function stops at quadruples.
#' Such an approach was proposed for these distinct starters in Colbourn (2004), who also provided
#' a smaller selection of starter vectors. (The approach is analogous to \code{\link{recBoseCA}}.)
#'
#' @references Colbourn (2004), Colbourn et al. (2006)
#'
#' @examples
#' ## mini example of Colbourn et al. (2006)
#' k <- 7; v <- 5
#' starter <- c(Inf, 0, 0, 1, 3, 2)
#' CS_CMMSSY(7, 5, starter=starter) ## 29 runs, 7 columns
#'
#' CS_CMMSSY(9, 6)  ## 46 rows
#'

#' @export
CS_CMMSSY <- function(k, v, start0=TRUE, starter=NULL, ...){
  Call <- sys.call()
  ## checks
  stopifnot(is.numeric(k), is.numeric(v))
  if (v == 2) stop(paste0("For v=2, use function KSK .\n It yields a strength 2 CA in ", N_KSK(k), "runs."))
  stopifnot(v >= 3)
  stopifnot(v %% 1 == 0)
  stopifnot(k > 2)
  stopifnot(is.logical(start0))
  if (!is.null(starter)){
    ## the starters must have one less element than the number of columns
    ## in this construction
    stopifnot(length(starter)>=k-1)
    stopifnot(starter[1]==Inf)
    stopifnot(all(starter[-1] %in% 0:(v-1)))
  }
  if (is.null(starter)){
    N <- N_CS_CMMSSY(t=2, k, v)
    if (is.na(N)) stop("This combination of k and v cannot be accommodated.")
    k_lookup <- k_CS_CMMSSY(N,v)  ## the column size corresponding to the N
    if (is.na(k_lookup)) stop("Unexpected error.")
    ## retrieve code that fetches starter, if possible
    starter=CMMSSYStarters[[as.character(v)]][[as.character(k_lookup-1)]]
    # zeile <- which(CMMSSYCombis$k==k_lookup &
    #                  CMMSSYCombis$v==v &
    #                  CMMSSYCombis$N==N)
    # stopifnot(length(zeile)==1)
    # ## conduct productPCA construction
    # code <- CMMSSYCombis$code[zeile]
    # #if (length(grep("starter", code, fixed=TRUE)) > 0)
    # return(eval(parse(text=CMMSSYCombis$code[zeile]))[,1:k])
  }

  ## create array for a specified starter
  ## and a built-in starter
  ## that was picked above

  G <- createG_CMMSSY(v)
  aus <- createCS_CMMSSY(starter, v, G)[,1:k]
  class(aus) <- c("ca", class(aus))
  aus <- CA_to_PCA(aus) ## for adding the attribute PCAstatus
  attr(aus, "Call") <- Call
  attr(aus, "t") <- 2
  attr(aus, "origin") <- "Cover starter CMMSSY 2006"
  if (!start0) aus <- aus + 1
  aus
}

#' @export
N_CS_CMMSSY <- function(t=2, k, v, theoretical=FALSE){
  ## checks
  stopifnot(is.numeric(k), is.numeric(v))
  stopifnot(k%%1==0, v%%1==0)
  stopifnot(is.logical(theoretical))
  if(k<2 ||  v<3 || !t==2) return(NA)
  if (theoretical) return((v-1)*k + 1)
  hilf <- CMMSSYCombis[CMMSSYCombis[,"v"]==v &
                                 CMMSSYCombis[,"k"]>=k ,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(min(hilf[,"N"]))
}

#' @export
k_CS_CMMSSY <- function(N, v){
  stopifnot(is.numeric(N), is.numeric(v))
  if (v == 2) return(NA)
  hilf <- CMMSSYCombis[CMMSSYCombis[,"v"]==v & CMMSSYCombis[,"N"]<=N,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(max(hilf[,"k"]))
}


createG_CMMSSY <- function(v){
  ## creates function that uses the additive group
  ## except for the Inf element
  ## levels are coded from 0 to v-2
  function(x) ifelse(x==Inf, Inf, (x+1) %% (v-1))
}

createCS_CMMSSY <- function(starter, v, G, ...){
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
  k <- length(starter) + 1

  ## first design row 1 (or 0) only
  ## then M
  ## then M with permuted elements (by G)
  mat <- M #  <- cbind(M, 1) ## SCA
  for (i in 1:(v-2)){
    M <- G(M)
    mat <- rbind(mat, M)
  }
  for (i in rev(c(0:(v-2),Inf))){
    mat <- rbind(i, mat)
  }
  ## last column
  mat <- cbind(mat, c(rep(0, v),
                      rep(1:(v-1), each=nrow(M))))
  mat[mat==Inf] <- v-1  ## fixe Inf values
  rownames(mat) <- NULL
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
