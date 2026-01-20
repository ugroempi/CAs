#' functions for cyclotomy constructions
#'
#' cyc provides cyclotomy constructions according to Colbourn (2010), which yield
#'  good high strength arrays for v=2,3,4,6 (see Details section); Ncyc and kcyc
#'  provide run sizes and maximum numbers of columns for the respective q, v and
#'  construction type.
#'
#' @rdname cyclotomy
#' @aliases cyclotomyCA
#' @aliases cyc
#' @aliases N_cyc
#' @aliases k_cyc
#'
#' @usage cyclotomyCA(t, k, v, ...)
#' @usage cyc(q, v, k=NULL, type=NULL, primitive=NULL)
#' @usage N_cyc(q, v, type)
#' @usage k_cyc(q, v, type)
#'
#' @param t requested strength (coverage)
#' @param k requested number of columns; for \code{cyc}, default \code{NULL} for maximum possible
#' @param v number of levels for each column
#' @param q prime or prime power with q=1(mod v)
#' @param type the construction type ("1","2","3","3a","3b","4a", or "4b"),
#' default implies "1"
#' @param primitive the primitive of the finite field to be used for the construction
#' (specification should not be needed, except for cases that are not adequately covered
#' by the internal data.frame \code{primedat})
#' @param ... currently not used
#'
#' @returns \code{cyclotomyCA} returns a matrix with \code{k} columns and values 0 to \code{v-1} in the
#' number of rows that can be obtained by the cyclotomy construction for strength \code{t};\cr
#' the matrix has class \code{ca}.\cr
#' @returns \code{cyc} returns a matrix with \code{k} columns and values 0 to \code{v-1} with the
#' number of runs corresponding to the construction (see Details section).
#' \code{k_cyc} and \code{N_cyc} return the maximum possible k and the run size N for the given parameters.
#'
#' @section Details:
#' Function \code{cyclotomyCA} allows to specify \code{t}, \code{k} and \code{v} and to
#' obtain the smallest \emph{implemented} CA from the cyclotomy construction
#' for this specification (if one exists); note that it is very well possible that there is
#' a smaller CA for the purpose from a cyclotomy construction. This package only automatically
#' implements those CAs that have optimum size as of the state of the Colbourn tables of best bounds
#' of November 2024; this is because various mistakes were found in tables of Colbourn (2010),
#' and the package author does not want to take responsibility for the strength of CAs according to
#' those tables from the paper. Users who know that a certain prime power \code{q}
#' yields a smaller CA for their requirements can use function \code{cyc} for creating that CA.\cr
#' The implemented constructions include three strength 5 cases attributed to
#' Torres-Jimenez, for which \code{q%%v==2} (and not 1, as required in Colbourn 2010).
#' Limited checks by sampling suggest that these are OK, but no full checks have been conducted,
#' and the implementation trusts that the Colbourn tables are correct (except in one case, which was
#' omitted because the alleged \code{q} is not a prime or prime power (13*229)).
#'
#' The function \code{cyc} is meant for experts and yields designs with columns in \code{v} levels based on
#' prime or prime power \code{q} as follows:
#'
#'   for N=k=q=1(mod v) prime power (Construction 4.1)\cr
#'   or N=q+v-1, k=q, with q=1(mod v) prime power (Construction 4.2)\cr
#'                (add constant rows with elements 1, 2, ..., v-1\cr
#'                 to construction 4.1, as these are missing)\cr
#'   or N=vq, k=q, with q=1(mod v) prime power (Construction 4.3)\cr
#'                (juxtapose vertically A, A+1, ..., A+v-1 (mod v);\cr
#'                 this is called developing modulo v)\cr
#'           or 4.3a: k=q+1\cr
#'                (extend A by column of 0s before developing)\cr
#'           or 4.3b: N=vq+v^2-v, k=q+1\cr
#'                (add non-zero constant rows (like in 4.2),\cr
#'                 then add column of 0s before developing)\cr
#'   or N=vq+v, k=q, with q=1(mod v) prime power (Construction 4.4)\cr
#'                add constant row of zeroes to A before developing\cr
#'           or 4.4a: k=q+1\cr
#'                also add constant column of zeroes\cr
#'           or 4.4b: N=vq+v^2-v, k=q+1 (same as 3.b)\cr
#'                add to A all nonzero constant rows\cr
#'                and a column of ones, then develop
#'
#' The parameters have to be specified by the user. They can, e.g.,
#' be selected from cyclotomy construction entries of \code{CAs:::colbournBigFrame},
#' or from the Colbourn (2010) paper.
#'
#' @section Warning: For valid parameters, function \code{cyc} always returns an array,
#' but coverage properties are only guaranteed for well-chosen parameters!
#'
#'
#' @examples
#' ## construction 1
#' ## 2-level
#' CA19.3.19.2 <- cyc(19, 2) ## q=19, k=19
#' coverage(CA19.3.19.2, 3)
#' coverage(CA19.3.19.2, 4)
#'
#' # 3-level
#' CA93.3.32.3 <- cyc(31, 3, type="3a")
#' coverage(CA93.3.32.3, 3)
#' dim(CA93.3.32.3)
#' k_cyc(31,3,"3a")
#' N_cyc(31,3,"3a")
#'
#' ## reduce k versus the default 32
#' CA93.3.23.3 <- cyc(31, 3, k=23, type="3a")
#' dim(CA93.3.23.3)
################################################################

#' @export
cyclotomyCA <- function(t, k, v, ...){
  Call <- sys.call()
  hilf <- CYCLOTOMYcat[CYCLOTOMYcat$t>=t & CYCLOTOMYcat$v==v & CYCLOTOMYcat$k>=k,]
  if (t > 6) stop("t is too large")
  if (v > 21) stop("v is too large")
  if (nrow(hilf)==0) stop("k is too large for the combination of t and v")
  N <- min(hilf$N)
  hilf <- hilf[which.min(hilf$N),,drop=FALSE]
  N <- hilf$N
  aus <- eval(parse(text=hilf$code))
  class(aus) <- c("ca", class(aus))
  attr(aus, "origin") <- paste0("Cyclotomy, q=", hilf$q, ", type=",
                                hilf$type)
  attr(aus, "Call") <- Call
  attr(aus, "t") <- t
  aus[,1:k]
}

#' @export
cyc <- function(q, v, k=NULL, type=NULL, primitive=NULL){
  ## k is the number of columns (default: as large as possible)
  ## v is the common number of levels
  ## q is the prime or prime power on which the construction is based
  ## type is the construction type ("1","2","3","3a","3b","4a","4b")

  if (is.null(type)) type <- "1"
  ## kcyc also checks its inputs
  hilf <- k_cyc(q, v, type)
  if (is.null(k)) k <- hilf else if (k>hilf) stop("k can be at most ", hilf)

  ## primedat is in sysdata.rda
  qs <- primedat$q
  firstprimitive <- primedat$primitive

  ## galois-field for prime
  gfsimple <- function(q){
    stopifnot(primes::is_prime(q))
    n <- 1
    p <- q
    xton <- 0
    inv <- NULL
    neg <- c(0, rev((q-1):1))
    root <- NULL
    plus <- outer(0:(q-1), 0:(q-1), function(X, Y) (X + Y) %% q)
    times <- outer(0:(q-1), 0:(q-1), function(X, Y) (X * Y) %% q)
    poly <- matrix(0:(q-1), q, 1)
    aus <- list(n=n, p=p, q=q, neg=neg, plus=plus, times=times, poly=poly)
    class(aus) <- c("Galoisfield", class(aus))
    aus
  }

  ## define minus for a gf
  gf_minus <- function(x,y,gf){
    ## calculates x - y, after reducing both by a mod operation to
    ## gf entries
    q <- gf$q; p <- gf$p
    x <- x%%q; y <- y%%q
    yn <- gf$neg[y+1]
    gf$plus[x+1,yn+1]
  }
  simple <- primes::is_prime(q)
  if (simple) gf <- gfsimple(q) else
  gf <- lhs::create_galois_field(q)
  ## start vector is the vector of "logarithms of 1:q" (mod v) w.r.t. the
  ##     base chosen as a primitive element (omega) of the group based on q
  ## retrieve primitive
  if (!is.null(primitive)) p <- primitive else p <- firstprimitive[which(qs==q)]
  ## create start vector

  xstart <- rep(0,q)
    cur <- p  ## current primitive
    for (i in 1:(q-1)) {
      xstart[cur+1] <- i  ## exponent of the primitive
      cur <- SOAs:::gf_prod(cur, p, gf)
    }
    if (!all(table(xstart)==1)) stop("wrong primitive")
    xstart <- xstart%%v
  ## creating the matrix A
  A <- matrix(NA,q,q)
  ## the complicated way of calculating j-i is needed for prime powers,
  ## otherwise (j-i)%%q would do
  ## might be worthwhile to separate, because the %% variant is faster
  if (!simple){
  for (i in 1:q)
    for (j in 1:q)
      A[i,j] <- xstart[gf_minus(j,i, gf)+1]
  }else{
    for (i in 1:q)
      for (j in 1:q)
        A[i,j] <- xstart[(j-i)%%q + 1]
  }
  if (type=="1") return(A[,1:k])
  if (type=="2") return(rbind(A,matrix(1:(v-1),nrow=v-1,ncol=q))[,1:k])
  if (type=="3") return(do.call(rbind, lapply(0:(v-1),
                        function(obj) (A+obj)%%v))[,1:k])
  if (type=="3a") return(do.call(rbind, lapply(0:(v-1),
                        function(obj) (cbind(A,0)+obj)%%v))[,1:k])
  if (type=="3b") return(do.call(rbind, lapply(0:(v-1),
                        function(obj) (cbind(rbind(A,
                                                   matrix(1:(v-1),
                                                          nrow=v-1,ncol=q)),
                                             0)+obj)%%v))[,1:k])
  if (type=="4") return(do.call(rbind, lapply(0:(v-1),
                        function(obj) (rbind(A,0)+obj)%%v))[,1:k])
  if (type=="4a") return(do.call(rbind, lapply(0:(v-1),
                        function(obj) (cbind(rbind(A,0),0)+obj)%%v))[,1:k])
  if (type=="4b") return(do.call(rbind, lapply(0:(v-1),
                        function(obj) (cbind(rbind(A,
                                                   matrix(1:(v-1),nrow=v-1,ncol=q)),
                                             1)+obj)%%v))[,1:k])
}

#' @export
## run size
N_cyc <- function(q, v, type){
  #if (!q%%v==1) stop("q mod v=1 is violated")
  if (!v<q) stop("v is too large")
  if (!type %in% c("1","2","3","3a","3b","4","4a","4b")) stop("invalid type")
  stopifnot(any(primes::is_prime(q^(1/(1:10)))))
  if (type=="1") return(q)
  if (type=="2") return(q+v-1)
  if (type %in% c("3","3a")) return(q*v)
  if (type %in% c("3b","4b")) return(v*(q+v-1))
  if (type %in% c("4", "4a")) return(v*(q+1))
}

#' @export
## maximum number of columns
k_cyc <- function(q, v, type){
  ## deactivated, because the Torres-Jimenez constructions
  ## for strength 5 have 2
  #if (!q%%v==1) stop("q mod v=1 is violated")
  if (!v<q) stop("v is too large")
  if (!type %in% c("1","2","3","3a","3b","4","4a","4b")) stop("invalid type")
  stopifnot(any(primes::is_prime(q^(1/(1:10)))))
  if (type %in% c("1","2","3","4")) return(q) else return(q+1)
}


