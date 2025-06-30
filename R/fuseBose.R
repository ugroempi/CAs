#' Function to create a strength 2 CA by fusion applied to a Bose OA
#'
#' A strength 2 CA in q+1 (q-1)-level columns and q^2-3 runs is created
#' from a q-level Bose array arranged as SCA
#'
#' @rdname fuseBose
#'
#' @aliases fuseBoseCA
#' @aliases fuseBose
#' @aliases N_fuseBose
#' @aliases k_fuseBose
#'
#' @usage fuseBoseCA(k=NULL, v=NULL, fixNA=TRUE, ...)
#' @usage fuseBose(q, fixNA=TRUE, ...)
#' @usage N_fuseBose(k=NULL, v=NULL, ...)
#' @usage k_fuseBose(N, v=NULL, ...)
#'
#' @param k \code{NULL} or requested number of columns (\code{q+1}); if \code{NULL}, inferred from \code{v} as \code{k=v+2}; if specified, \code{k}-1 must be a prime power, or \code{k} will be increased to the nearest such number.
#' @param v \code{NULL} or number of levels; if \code{NULL}, inferred from \code{k} as \code{v=k+2} or from \code{N} as \code{sqrt(N+3)-1}; if \code{k} is also specified, it takes precedence.
#' @param q a prime power
#' @param fixNA logical; if FALSE, keeps flexible values; if TRUE, the first \code{q}-3 rows are made constant.
#' @param k \code{NULL} or requested number of columns (\code{q+1}); if \code{NULL}, inferred from \code{v} as \code{k=v+2}; if specified, \code{k}-1 must be a prime power, or \code{k} will be increased to the nearest such number.
#' @param v \code{NULL} or number of levels; if \code{NULL}, inferred from \code{k} as \code{v=k+2} or from \code{N} as \code{sqrt(N+3)-1}; if \code{k} is also specified, it takes precedence.
#' @param N maximum acceptable number of runs
#' @param ... currently not used
#'
#' @details
#' Colbourn (2008) proposed this particular variant of fusion, which allows to
#' reduce the number of runs by 3 for strength 2 and at most q+1 columns.
#' It creates a CA(q^2-3, 2, q+1, q-1) from an OA(q^2, 2, q+1, q).
#' The reduction in rows is achieved by using a graph coloring approach
#' (Lemma 3.2 of Colbourn 2008).
#'
#' @returns \code{fuseBose} returns a (\code{q^2-3} x \code{q+1} matrix with columns in \code{q}-1 levels labeled from 0 to q-2 that is a covering array of strength 2.\cr
#' \code{N_fuseBose} and \code{k_fuseBose} return a named vector with elements \code{N}, \code{k}, \code{v} and \code{q}.
#'
#' @references Colbourn (2008, Lemma 3.2)
#'
#' @examples
#' # create a CA(22,2,6,4)
#' fuseBose(5)  ## two constant rows
#' fuseBose(5, fixNA=FALSE) ## two flexible values
#' head(fuseBoseCA(6, 4)) ## two constant rows
#' attributes(fuseBoseCA(6, 4, fixNA=FALSE)) ## flexible values
#'
#' # create a CA(46,2,8,6)
#' dim(D <- fuseBose(7))
#' coverage(D, 2)
#' eCAN(2,8,6) ## four runs difference
#'
#' # create a CA(61,2,9,7)
#' dim(fuseBose(8))
#' eCAN(2, 9, 7)  ## two runs difference
#'
#' # create a CA(166,2,14,12)
#' dim(fuseBose(13))
#' eCAN(2, 14, 12)  ## one run difference
#'
#' # querying N
#' N_fuseBose(k=21)
#' eCAN(2, 24, 22) ## one run difference
#' eCAK(2, 526, 22) ## as good as it gets
#'
#' N_fuseBose(v=6)
#'

#' @export
fuseBoseCA <- function(k=NULL, v=NULL, fixNA=TRUE, ...){
  Call <- sys.call()
  hilf <- N_fuseBose(k,v,...)
  aus <- fuseBose(hilf["q"], fixNA=fixNA, ...)
  class(aus) <- c("ca", class(aus))
  attr(aus, "Call") <- Call
  attr(aus, "origin") <- "Colbourn 2008, Lemma 3.2"
  if (!fixNA){
    if (any(is.na(aus)))
      attr(aus, "flexible") <- list(value=NA, profile=colSums(is.na(aus)))
  }
  aus
}

  #' @export
fuseBose <- function(q, fixNA=TRUE, ...){
  stopifnot(is.logical(fixNA))
  ## first row immediately dropped
  arr <- SCA_Bose(q)[-1,]
  arr[arr==0] <- NA ## removing the level 0 entirely
  ## arrsub has the first rows constant in increasing order
  arrsub <- arr[,-(q+1)]

  vertices <- which(is.na(arrsub), arr.ind=TRUE)
  ## store the NA-positions, which are graph vertices

  G <- igraph::make_empty_graph(n=nrow(vertices), directed = FALSE)

  ## find edges
  paare <- nchoosek(q, 2)
  for (paar in 1:ncol(paare)){
    spalten <- paare[,paar]
    for (i in 1:2){
      ## loop over the to-be-deleted rows
      ## assign edges as indicated in Colbourn 2008 (proof of Lemma 3.2)
      row1 <- which(arrsub[,spalten[1]]==i & is.na(arrsub[,spalten[2]]))
      row2 <- which(is.na(arrsub[,spalten[1]]) & arrsub[,spalten[2]]==i)
      from <- which(vertices[,1]==row1 & vertices[,2]==spalten[2])
      to <- which(vertices[,1]==row2 & vertices[,2]==spalten[1])
      G <- igraph::add.edges(G, c(from, to))
    }
  }
  ## yields numbers 1 and 2, perfect for being assigned to
  ## NA values
  farben <- igraph::greedy_vertex_coloring(G)
  ## suitably assign values in order to achieve coverage
  for (i in 1:nrow(vertices)){
    arrsub[vertices[i,1], vertices[i,2]] <- farben[i]
  }
  aus <- cbind(arrsub[-(1:2),], arr[-(1:2),q+1])-1
  if (fixNA) aus[1:(q-3), q+1] <- aus[1:(q-3), 1]
  aus
}

#' @export
N_fuseBose <- function(k=NULL, v=NULL, ...){
  primpotenzen <- primedat$q
  if (is.null(k) && is.null(v)) stop("k and v must not both be NULL")
  if (!is.null(k)){
    stopifnot(is.numeric(k))
    stopifnot(k %% 1==0)
    if (!is.null(v)){
      stopifnot(k <= v+2)
      ## check for prime
      if (!(v+1) %in% primedat$q) return(NA)
      kalt <- k
      q <- v+1
    }
    else{
      v <- k-2
      kalt <- k
      if (!(v+1) %in% primedat$q){
        # increase k until v+1 is prime power
        while(!(k-1) %in% primedat$q){
          k <- k+1
          v <- k-2
        }
      }
      q <- v+1
    }
  }else{
    ## now k is NULL
    stopifnot(is.numeric(v))
    stopifnot(v %% 1==0)
    ## check for prime
    if (!(v+1) %in% primedat$q) return(NA)
      # stop("v+1 must be a prime power")
    q <- v+1
    k <- kalt <- q+1   ## k was NULL
  }
  c(N=q^2-3, k=kalt, kmax=q+1, v=v, q=q)
}

### #' do not @export yet
k_fuseBose <- function(N, v=NULL, ...){
  primpotenzen <- primedat$q
  hilf <- primpotenzen^2 - 3
  if (!(N %in% (hilf))){
    N <- max(hilf[hilf <= N ])
    message("Invalid N was reduced to closest valid N")
  }
  q <- as.integer(sqrt(N+3))  ## must be prime power
  if (!is.null(v)){
    stopifnot(is.numeric(v))
    stopifnot(v %% 1==0)
    if (!(v+1) %in% primedat$q)
      message("v is invalid and was ignored.")
    if (!(v+1 == q))
      message("v is not compatible with N and was ignored.")
    }
   v <- q-1
   c(N=q^2-3, k=q+1, v=v, q=q)
}

