#' Function to create a strength 2 CA by projection of a Bose OA
#'
#' A strength 2 CA with q^2-c runs with q+1 (q-c)-level columns and c s-level columns is created from a q-level Bose array arranged as SCA
#'
#' @rdname projectionBose
#'
#' @aliases projBoseCA
#' @aliases projectionBose
#' @aliases N_projectionBose
#' @aliases k_projectionBose
#'
#' @usage projBoseCA(k, v, fixNA=TRUE, ...)
#' @usage projectionBose(q, c=1, s=q-c, fixNA=FALSE, ...)
#' @usage N_projectionBose(k, v, cmax=Inf, ...)
#' @usage k_projectionBose(N, v, cmax=Inf, ...)
#'
#' @param k requested number of columns
#' @param v number of levels
#' @param q a prime power
#' @param c an integer number, 0 <= c <= \code{q-2}; for \code{c=0}, the unchanged array created by \code{\link{SCA_Bose}} is returned.
#' @param s an integer number, 2 <= s <= \code{q-c}; for \code{c=0}, \code{s} is irrelevant.
#' @param fixNA logical; if FALSE, keeps flexible values; if TRUE, flexible values are arbitrarily fixed.
#' @param N maximum acceptable number of runs
#' @param cmax the maximum number permitted for \code{c}; setting \code{cmax} prevents ridiculous values of \code{N} for cases
#' for which \code{k} and \code{v} are too far apart
#' @param ... currently not used
#'
#' @returns \code{projBose} returns a covering array in \code{k} columns at \code{v} levels each with \code{v} constant rows,
#' if this is reasonably possible by projection of a strength 2 Bose OA, which should be the case whenever
#' \code{v<k} and \code{k-v} not larger than seven and \code{v+(k-v-1)/2} a prime power.\cr
#' @returns \code{projectionBose} returns a (\code{q^2-c} x \code{q+1+c} matrix with the first \code{q+1} columns
#' in \code{q-c} levels (labeled from 0 to \code{q-c-1}) and the last \code{c} columns in \code{s} levels (labeled from 0 to \code{s-1})
#' that is a covering array of strength 2.\cr
#' \code{N_projectionBose} and \code{k_projectionBose} return a named
#' vector with elements \code{N}, \code{k}, \code{kmax}, \code{v}, \code{q} and \code{c}.
#'
#' @references Colbourn (2008, Theorem 2.3), Torres-Jimenez et al. (2019, Fig. 25)
#'
#' @examples
#' # create a CA(24,2,7,4)
#' projectionBose(5)  ## four constant rows, some flexible values
#' projectionBose(5, fixNA=TRUE) ## flexible values arbitrarily filled
#' projectionBose(5, c=2, s=2) ## 6 3-level columns (3=5-2), two 2-level columns
#'
#' # create a CA(48,2,9,6)
#' dim(D <- projectionBose(7))
#' coverage(D, 2)
#' eCAN(2,9,6) ## two runs difference
#'
#' # create a CA(63,2,9,7)
#' dim(projBoseCA(9,7))
#' eCAN(2, 9, 7)  ## four runs difference
#'
#' # create a CA(166,2,15,12)
#' dim(projBoseCA(15,12))
#' eCAN(2, 15, 12)  ## optimal
#'
#' # querying k or N
#' N_projectionBose(k=21, v=15)
#' k_projectionBose(N=357, v=15)
#' eCAN(2, 21, 15) ## the above design is the basis
#' eCAK(2, 357, 15) ## four more columns
#'

#' @export
projBoseCA <- function(k, v, fixNA=TRUE, ...){
  Call <- sys.call()
  stopifnot(v < k)
  hilf <- N_projectionBose(k, v)  ## with default cmax=Inf
  if (is.na(hilf[1])) stop("This construction is not available for this combination of k and v.")
  aus <- projectionBose(hilf["q"], hilf["c"], fixNA=fixNA, ...)[,1:k]
  class(aus) <- c("ca", class(aus))
  attr(aus, "PCAstatus") <- list
  attr(aus, "call") <- Call
  attr(aus, "origin") <- "projection (Colbourn 2008 Thm 2.3)"
  aus
}

#' @export
projectionBose <- function(q, c=1, s=q-c, fixNA=FALSE, ...){
  ## s is the number of levels of the additional columns
  ## if it equals q-c, the resulting CA is uniform
  stopifnot(is.numeric(q), is.numeric(c), is.numeric(s))
  stopifnot(c>=0)
  stopifnot(c<=q-2)    ## make sure at least two levels remain
  stopifnot(s <= q-c)  ## make sure that s is permissible
  stopifnot((c%%1)==0)
  stopifnot(is.logical(fixNA))
  arr <- SCA_Bose(q)                  ## also checks for appropriate q
  if (c==0) return(arr)
  newmat <- cbind(arr, matrix(NA, nrow(arr), c))
       ## c added columns, initialized as flex

  ## immediately use final symbols only for the first s columns
  replacesel <- setdiff(0:(q-1), 0:(c-1))
  for (j in 1:c){
    ## loop through the first c rows (rows to be removed)

    for (k in 1:s){
      ## determine C_k
      hilf <- setdiff(which(newmat[,k]==newmat[j,k]), (1:c))
      ## fill the new column
      newmat[hilf, q+1+j] <- k - 1 + c  ## immediately correct levels
                                        ## no after-treatment needed
                                        ## except subtracting c
      ## treat the first s columns
      if (length(replacesel) < length(hilf))
        replace <- sample(
          c( replacesel, rep(NA, length(hilf)-length(replacesel)) )
          )
      else replace <- sample(replacesel)
             ## permutation of target values (starting with c)
      newmat[hilf, k] <- replace
    }

    ## treat value j-1 in columns s+1 to q+1 (wild card)
    newmat[(c+1):(q^2),(s+1):(q+1)][which(newmat[(c+1):(q^2),(s+1):(q+1)]==j-1)] <- NA
  }

  aus <- newmat[-(1:c),] - c
  ## make first rows constant
  for (i in 1:(q-c)){
     if (s==q-c || i <= s)
       aus[i, (q+1):(q+1+c)] <- i-1
     else
       aus[i, (q+1)] <- i-1
  }
  missings1 <- which(is.na(aus[,1:(q+1)]))
  missings2 <- which(is.na(aus[,(q+2):(q+1+c)]))
  if (fixNA){
    aus[,1:(q+1)][missings1] <- sample(0:(q-c-1), length(missings1), replace = TRUE)
    aus[,(q+2):(q+1+c)][missings2] <- sample(0:(s-1), length(missings2), replace = TRUE)
  }
  aus
}

#' @export
N_projectionBose <- function(k, v, cmax=Inf, ...){
  ## k and v must be given
  ## k is a minimum, v must be matched exactly
  ## there must be a prime power q with v+c=q and k<=q+1+c
  primpotenzen <- primedat$q
    stopifnot(is.numeric(k), is.numeric(v), is.numeric(cmax))
    stopifnot(k %% 1==0, v%%1==0)
    kalt <- k
    if (v > k) k <- v+1  ## prevent nonsensical results
    if (!((k-v)%%2)==1){
      k <- k+1
      }
    c <- (k-v-1)/2
    q <- NA
    ## allows q=v+c>k, since there are instances where this is the current best array
    ## in the Colbourn tables
    while(v+c <= k + 1){
    if ((v+c) %in% primpotenzen){
      q <- v+c
      break
      }else
      c <- c+1
      }
    if (is.na(q)) return(NA) ## no adequate scenario found
    if (c > cmax) return(NA) ## no adequate and permissible scenario found
  c(N=q^2-c, k=kalt, kmax=q+1+c, v=v, q=q, c=c)
}

#' @export
k_projectionBose <- function(N, v=NULL, cmax=Inf, ...){
  ## v=NULL is interpreted as 2
  primpotenzen <- primedat$q
  if (is.null(v)) v <- 2
  stopifnot(is.numeric(v), is.numeric(N))
  stopifnot(N%%1==0, v%%1==0)
  hilfc <- primpotenzen - v
  hilf <- primpotenzen^2 - hilfc

  if (!(N %in% hilf)){
    pos <- which.max(hilf[hilf <= N ])
    N <- hilf[pos]
    c <- hilfc[pos]
    message("Invalid N was reduced to closest valid N")
  }else{
    pos <- which(hilf==N)
    c <- hilfc[pos]
  }
  if (c<0 || c>cmax) return(NA)
  q <- v + c
  c(N=q^2-c, k=q+1+c, v=v, q=q, c=c)
}

