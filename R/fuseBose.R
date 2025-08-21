#' Function to create a strength 2 CA by fusion applied to a Bose OA
#'
#' A strength 2 CA in q+1 (q-1)-level columns and q^2-3 runs is created
#' from a q-level Bose array arranged as SCA
#'
#' @rdname fuseBose
#'
#' @aliases fuseBoseCA
#' @aliases fuseBose
#' @aliases fuseBushtCA
#' @aliases N_fuseBose
#' @aliases k_fuseBose
#' @aliases N_fuseBusht
#' @aliases k_fuseBusht
#'
#' @usage fuseBoseCA(k=NULL, v=NULL, fixNA=TRUE, ...)
#' @usage fuseBose(q, fixNA=TRUE, ...)
#' @usage fuseBushtCA(t, k, v, fixNA=TRUE, ...)
#' @usage N_fuseBose(k=NULL, v=NULL, ...)
#' @usage k_fuseBose(N, v=NULL, ...)
#' @usage N_fuseBusht(t, k, v, ...)
#' @usage k_fuseBusht(t, N, v, ...)
#'
#' @param k \code{NULL} or requested number of columns (\code{q+1}); if \code{NULL}, inferred from \code{v} as \code{k=v+2}; if specified, \code{k}-1 must be a prime power, or \code{k} will be increased to the nearest such number.
#' @param v \code{NULL} or number of levels; if \code{NULL}, inferred from \code{k} as \code{v=k+2} or from \code{N} as \code{sqrt(N+3)-1}; if \code{k} is also specified, it takes precedence.
#' @param q a prime power
#' @param t the requested strength
#' @param fixNA logical; if FALSE, keeps flexible values; if TRUE, the first \code{q}-3 rows are made constant.
#' @param k \code{NULL} or requested number of columns (\code{q+1}); if \code{NULL}, inferred from \code{v} as \code{k=v+2}; if specified, \code{k}-1 must be a prime power, or \code{k} will be increased to the nearest such number.
#' @param v \code{NULL} or number of levels; if \code{NULL}, inferred from \code{k} as \code{v=k+2} or from \code{N} as \code{sqrt(N+3)-1}; if \code{k} is also specified, it takes precedence.
#' @param N maximum acceptable number of runs
#' @param ... currently not used
#'
#' @details
#' Colbourn (2008) proposed the particular variant of fusion for a Bose CA,
#' which allows to reduce the number of runs by 3 for strength 2 and at most q+1 columns.
#' It creates a CA(q^2-3, 2, q+1, q-1) from an OA(q^2, 2, q+1, q).
#' The reduction in rows is achieved by using a graph coloring approach
#' (Lemma 3.2 of Colbourn 2008).
#'
#' The fusion of Bush strength \code{t} CAs is included for cases
#' for which no other better CA is implemented. It creates a
#' CA(q^t-2*c, t, q+1, q-c) (with c the distance of v to the next largest prime power q)
#'
#' @returns \code{fuseBose} returns a (\code{q^2-3} x \code{q+1} matrix with columns in \code{q}-1 levels labeled from 0 to q-2 that is a covering array of strength 2.\cr
#' \code{fuseBoseCA} returns a matrix of class \code{ca} from inputs \code{k} and \code{v}.\cr
#' \code{fuseBushtCA} returns a matrix of class \code{ca} which is obtained from the smallest possible strength \code{t} Bush construction SCA
#' created with \code{\link{SCA_Busht}} by the necessary number of fuse steps.\cr
#' \code{N_fuseBose} and \code{k_fuseBose}, \code{N_fuseBusht} and \code{k_fuseBusht}
#' return a named vector with elements \code{N}, \code{k}, \code{v} and \code{q}.
#'
#' @references Colbourn (2008, Lemma 3.2) for \code{fuseBose}
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
#' # create a CA(508,3,9,6)
#' dim(fuseBushtCA(3, 9, 6))
#' eCAN(2, 14, 12)  ## one run difference
#'
#' # querying N
#' N_fuseBose(k=21)
#' eCAN(2, 24, 22) ## one run difference
#' eCAK(2, 526, 22) ## as good as it gets
#'
#' N_fuseBose(v=6)
#'
#' ## run size from fuseBusht for t=3,k=9,v=7
#' N_fuseBusht(3,9,7)
#'
#' ## fuseBusht is very close to eCAN
#' ## postopNCK removes a few rows
#' ## (not in its implementation in this package)
#' Ns(3, 18, 15)
#'
#' ## fuseBusht is the best implemented solution
#' Ns(3,14,12)
#' ## it is behind the eCAN entry
#' eCAN(3, 14, 12)
#' ## apparently, postopNCK can remove 25 runs
#' ## with the implementation in the package the effort
#' ## for this appears somewhat unreasonable

#' @export
fuseBoseCA <- function(k=NULL, v=NULL, fixNA=TRUE, ...){
  Call <- sys.call()
  hilf <- N_fuseBose(k,v,...)
  aus <- fuseBose(hilf["q"], fixNA=fixNA, ...)[,1:k]
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
fuseBushtCA <- function(t, k, v, fixNA=TRUE, ...){
  Call <- sys.call()
  hilf <- try(N_fuseBusht(t,k,v,...))
  if ("try-error" %in% class(hilf))
    stop("fuseBushtCA does not work for this setting")
  if (hilf["q"]==hilf["v"]){
    aus <- SCA_Busht(v, t)
    attrs <- attributes(aus)
    if (k < ncol(aus)){
        aus <- aus[,1:k]
        attrs$dim[2] <- k
        attrs$PCAstatus$k2 <- max(0, k-attrs$PCAstatus$k1)
        attrs$PCAstatus$k1 <- min(k, attrs$PCAstatus$k1)
        attributes(aus) <- attrs
    }
    return(aus)
  }
  ## now fusion takes place
  q <- hilf["q"]
  aus <- SCA_Busht(q, t)
  attrs <- attributes(aus)
  attrs$PCAstatus <- NULL
  attrs$dim[1] <- attrs$dim[1] - 2*(q-v)
  attrs$origin <- c(attrs$origin, paste(rep("fuse", q-v), collapse=" "))
  aus <- fuse(aus, vprior=q, vpost=v, fixNA=fixNA, ...)[,1:k]
  attrs$dim <- dim(aus)
  attributes(aus) <- attrs
  attr(aus, "Call") <- Call  ## overwrite
  if (!fixNA){
    if (any(is.na(aus)))
      attr(aus, "flexible") <- list(value=NA, profile=colSums(is.na(aus)))
  }
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

#' @export
k_fuseBose <- function(N, v=NULL, ...){
  ## 3 runs can be omitted
  q <- NA
  stopifnot(is.numeric(N))
  stopifnot(N%%1==0)
  if (!is.null(v)){
    stopifnot(is.numeric(v))
    stopifnot(v %% 1==0)
    if (!(v+1) %in% primedat$q)
      stop("v+1 must be a prime or prime power")
    else q <- v+1
  }
  if (!is.na(q)){
    if (N < q^2 - 3) return(NA)
    return(c(N=q^2-3, k=q+1, v=v, q=q))
  }
  ## now unspecified v
  hilf <- primedat$q^2 - 3
  if (!(N %in% (hilf))) N <- max(hilf[hilf <= N ])
  q <- as.integer(sqrt(N+3))  ## must be prime power
   v <- q-1
   c(N=q^2-3, k=q+1, v=v, q=q)
}

#' @export
N_fuseBusht <- function(t, k, v, ...){
  #if (t==2) return(NA)
  primpotenzen <- primedat$q
  stopifnot(is.numeric(k))
  stopifnot(k %% 1==0)
  ## find next largest suitable prime
  q <- min(primpotenzen[which(primpotenzen >= max(v, k-2))])
  if (q==k-2){
    ## this is suitable for powers of 2 with t=3 only
    if (!(((q/2)%%1)==0 && t==3)) {
      primpotenzen <- setdiff(primpotenzen, q)
      q <- min(primpotenzen[which(primpotenzen >= max(v, k-2))])
    }
  }
  ## now q is permissible
  kalt <- k
  ## kmax depends on the q,t combination
  kmax <- ifelse(t==3 && ((q/2)%%1)==0, q + 2, q + 1)

  ## prevent cases for which fewer than 2 levels would remain
  stopifnot(q-v <= v-2)
  c(N=q^t-2*(q-v), k=kalt, kmax=kmax, v=v, q=q)
}

#' @export
k_fuseBusht <- function(t, N, v, ...){
  primpotenzen <- primedat$q
  stopifnot(is.numeric(N))
  stopifnot(N %% 1==0)
  ## find next largest suitable prime
  q <- min(primpotenzen[which(primpotenzen >= max(v, N^(1/t)))])
  Nalt <- N
  c(N=q^t-2*(q-v), k = k_SCA_Busht(t, q^t, q), Nmax=Nalt, v=v, q=q)
}

