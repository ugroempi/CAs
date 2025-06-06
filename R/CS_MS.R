### Construction of Meagher and Stevens 2005 based on subgroup of symmetric group with fixed 0
### and circular matrix created from a cover starter vector
###    taken fom MeagherStevensStarters
#' Strength 2 CAs constructed by cover starters (Meagher and Stevens)
#'
#' Function to construct a strength 2 CA based on cover starters according to Meagher and Stevens.
#' k and v must be specified, the number of runs is ((v-1)*k+1) (possibly for a larger value of k
#' if the construction does not exist for the requested k)
#'
#' @rdname CS_MS
#'
#' @aliases CS_MS
#' @aliases SCA_MS
#' @aliases N_CS_MS
#' @aliases k_CS_MS
#'
#' @usage CS_MS(k, v, start0=TRUE, starter=NULL, ...)
#' @usage SCA_MS(k, v, start0=TRUE, starter=NULL, ...)
#' @usage N_CS_MS(t=2, k, v, theoretical=FALSE)
#' @usage k_CS_MS(t=2, N, v)
#'
#' @param k number of factors
#' @param v number of levels, \code{v} >= 3
#' @param start0 should the output design's levels start with zero (default) or 1?
#' @param starter an optional cover starter vector.\cr
#' This may be useful for cases with k larger than implemented in the stored cover starters.\cr
#' CAUTION: Responsibility for achieving the coverage rests with the user, please check!
#' @param t requested strength (must be 2, added of unified interface with other N_ and k_ functions)
#' @param theoretical logical: default requests size available from using stored cover starters,\cr
#' change to \code{TRUE} for the size assuming existence of a cover starter for \code{k} and \code{v}
#' @param N affordable number of runs
#' @param ... currently not used
#'
#' @returns \code{CS_MS} returns the smallest possible strength 2 CA from the Meagher/Stevens construction for \code{k} experimental factors in \code{v} levels; it is an \code{N x k} matrix with entries from 0 to \code{v-1} or from 1 to \code{v}, depending on \code{start0}.\cr
#' \code{SCA_MS} does the same, but arranges the result as an SCA (see Section "Partitioned Covering Array (PCA)" in \code{\link{productPCA}}).\cr
#' \code{N_CS_MS} returns the minimum run size needed for the Meagher/Stevens construction with \code{k} factors at \code{v} levels, or \code{NA}, if no such construction exists.\cr
#' \code{k_CS_MS} returns the maximum number of \code{v}-level factors that can be accommodated in at most \code{N} runs, or \code{NA}, if no such construction exists.
#'
#' @section Details:
#' Meagher and Stevens (2005) denote designs with the factors as the rows and the runs as the columns, this package does the opposite.
#'
#' The construction uses a starter vector of \code{k} elements with first element zero and subsequent elements from 1 to \code{v-1}.
#' These were identified by searches (exhaustive for \code{v <= 10} with \code{v+1 <= k <= 2*(v+1)},
#' hill-climbing for \code{11<=k<=20}). Starter vectors identified by Meagher (2005b) in her thesis,
#' and in part only provided in the unpublished technical report Meagher (2005a),
#' are collected in the internal list \code{MeagherStevensStarters},
#' whose list entry for each \code{v} is a list of starter vectors named with \code{k};
#' search-based vectors for small \code{v} (3, 4, 5) have been added.\cr
#' Note that there can be larger, and in some cases also smaller, starter vectors than the ones listed, which is why function \code{CS_MS}\cr
#' permits the specification of a starter vector. \cr
#' In search of a starter vector for \code{k+1} if a starter vector for \code{k} is available,
#' inserting a "1" in the second position is often successful (this technique has also been used for finding suitable vectors for cases for which Meagher only asserted existence without providing a vector).
#'
#' The starter vector is cycled into a \code{k x k} matrix \code{M} (internal function \code{circular}),
#' which is one key ingredient of the construction. The other key ingredient is the subgroup \code{G} of
#' the symmetric group of permutations of \code{v} elements, which cyclically permutes the second to last of the \code{v} elements and leaves the first element fixed (i.e., has \code{v-1} elements).
#' The resulting matrix consists in horizontally concatenating a vector of zeroes with the
#' \code{v-1} results of applying to \code{M} the functions corresponding to a group element
#' of \code{G}.\cr
#' Group manipulations rely on R package \pkg{permutations} (which uses coding 1 to \code{v}).
#'
#' If there is no starter vector for \code{k} in the list for \code{v}, the function tries to use the next largest \code{k}, if possible. If that is not possible, an error is thrown.
#'
#' Functions \code{N_CS_MS} and \code{k_CS_MS} return the minimum \code{N} or the maximum
#' \code{k} for the requested scenario. In case of \code{N_CS_MS}, the default to use the
#' available starters can be overridden by the argument \code{theoretical}, which allows
#' to inspect how many runs can be obtained by a suitable user-specified starter vector.
#'
#' @references Meagher and Stevens (2005), Meagher (2005a), Meagher (2005b)
#'
#' @examples
#' ## mini example of Meagher / Stevens
#' k <- 5; v <- 3
#' # starter <- c(0,1,1,1,2)
#' CS_MS(5, 3) ## 11 runs
#'
#' ## smallest element of Meagher / Stevens table II
#' ## k = 9, v = 6
#' D <- CS_MS(9, 6)
#' dim(D)
#' coverage(D,2)
#' ## the matrix M from cyclic permutation of the starter vector
#' CAs:::MeagherStevensStarters[["6"]][["9"]]
#' ## is treated by cyclic permutation of all levels expect the lowest
#' ## (v-1) k x k matrices
#' ## and a constant column with the lowest level is prepended
#'
#' ## use k=8, which is not directly in the starter table
#' D <- CS_MS(8, 6)
#' dim(D)    ## works by finding the smallest k that is in the table for v=6
#' N_CS_MS(k=8, v=6)
#' N_CS_MS(k=8, v=6, theoretical=TRUE)  ## in this case unachievability is proven
#' N_CS_MS(k=15, v=6) ## not implemented with stored cover starter
#' N_CS_MS(k=15, v=6, theoretical=TRUE) ## if a cover starter is found and provided,
#'                                  ## 76 runs are possible with the construction
#'
#' try(D <- CS_MS(14, 6))  ## not implemented
#'                         ## (though a starter might exist,
#'                         ## but likely unattractive result)
#' N_CS_MS(k=15, v=6, theoretical=TRUE) ## 76 runs possible, if suitable starter provided
#' CAs:::MeagherStevensStarters[["6"]][["14"]] ## trying a guess for a starter
#' starter15 <- c(0,1,1,1,1,1,1,1,1,1,2,4,3,1,2)
#' D <- CS_MS(15, 6, starter=starter15)
#' coverage(D, 2)          ## successful guess
#' dim(D)

#'@export
CS_MS <- function(k, v, start0=TRUE, starter=NULL, ...){
  ## checks
  stopifnot(is.numeric(k), is.numeric(v))
  if (v == 2) stop(paste0("For v=2, use function KSK .\n It yields a strength 2 CA in ", N_KSK(k), "runs."))
  stopifnot(v >= 3)
  stopifnot(v %% 1 == 0)
  stopifnot(k > 2)
  stopifnot(is.logical(start0))
  if (!is.null(starter)){
    stopifnot(length(starter)==k)
    stopifnot(starter[1]==0)
    stopifnot(all(starter[-1] %in% 1:(v-1)))
  }
  if (is.null(starter)){
    N <- N_CS_MS(t=2, k, v)
    if (is.na(N)) stop("This combination of k and v cannot be accommodated.")
    k_lookup <- k_CS_MS(t=2,N,v)
    if (is.na(k_lookup)) stop("Unexpected error.")
    ## retrieve starter
    vc <- as.character(v)
    kc <- as.character(k_lookup)
    starter <- MeagherStevensStarters[[vc]][[kc]]
  }
  ## create array
  aus <- createCS_MS(starter, v, createG(v))[,1:k]
  if (start0) aus <- aus - 1
  aus
}

#'@export
N_CS_MS <- function(t=2, k, v, theoretical=FALSE){
  if (!t==2) return(NA)
  if (!v>=3) return(NA)
  ## checks
  stopifnot(is.numeric(k), is.numeric(v))
  stopifnot(k>=2, v>=3)
  stopifnot(k%%1==0, v%%1==0)
  stopifnot(is.logical(theoretical))
  if (v == 2) stop(paste0("For v=2, use function KSK .\n It yields the optimum run size, which is ", N_KSK(k), "."))
  if (k <= v+1) message("For k <= v+1, Meagher/Stevens may not be the best construction.")
  if (theoretical) return((v-1)*k + 1)
  hilf <- MeagherStevensCombis[MeagherStevensCombis[,"v"]==v &
                                 MeagherStevensCombis[,"k"]>=k ,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(min(hilf[,"N"]))
}

#'@export
k_CS_MS <- function(t=2, N, v){
  if (!t==2) return(NA)
  stopifnot(is.numeric(N), is.numeric(v))
  if (v == 2) stop(paste0("For v=2, use function KSK .\n It yields the maximum number of factors, which is ", k_KSK(N), "."))
  hilf <- MeagherStevensCombis[MeagherStevensCombis[,"v"]==v & MeagherStevensCombis[,"N"]<=N,,drop=FALSE]
  k_cand <- floor((N-1)/(v-1))
  if (k_cand > 2*(v+1)) k_cand <- 2*(v+1)
  if (nrow(hilf)==0) return(NA)
  else return(max(hilf[,"k"]))
}

### auxiliary functions

circular <- function(vec){
  ## matrix of cyclic permutations
  ## used in createCS_MS
  k <- length(vec)
  aus <- matrix(NA,k,k)
  hilf <- vec
  aus[,1] <- hilf
  for (i in 2:k){
    hilf <- c(hilf[k], hilf[1:(k-1)])
    aus[,i] <- hilf
  }
  aus
}

#'@importFrom permutations as.cycle
#'@importFrom permutations as.function.permutation
createG <- function(v){
  ## creates function that permutes levels 2:v and keeps level 1 fixed
  ## requires levels coded as 1 to v
  wordmat <- rbind(c(1,3:v,2))
  Gperm <- permutations::as.cycle(wordmat)
  ## make function for group action
  permutations::as.function.permutation(Gperm)
}

createCS_MS <- function(starter, v, G, ...){
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
  M <- circular(starter) + 1
  k <- length(starter)

  ## first design row 1 (or 0) only
  ## then M
  ## then M with permuted elements (by G)
  ## create in transposed form (as used in the Meagher Stevens paper)
  mat <- cbind(1,M)
  for (i in 2:(v-1)){
    M <- matrix(G(M), nrow=k)
    mat <- cbind(mat, M)
  }
  ## return with columns for factors and rows for runs
  t(mat)
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

#' @export
SCA_MS <- function(k, v, start0=TRUE, starter=NULL, ...){
  hilf <- CS_MS(k, v, start0=start0, starter=starter, ...)
  CA_to_PCA(hilf)
}
