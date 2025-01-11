#' function for cyclotomy constructions
#'
#' provides cyclotomy constructions according to Colbourn (2010), which yield
#'  good high strength arrays for v=2,3,4,6 (see Details section)
#'
#' @param k number of columns
#' @param t requested strength (\code{t} < \code{k})
#' @param v number of levels for each column
#' @param q prime or prime power with q=1(mod v)
#' @param type the construction type ("1","2","3","3a","3b","4a", or "4b"),
#' default implies "1"
#' @param primitive the primitive of the finite field to be used for the construction
#' (specification should not be needed, except for cases that are not adequately covered
#' by the internal data.frame \code{primedat})
#'
#' @return a matrix with \code{k} columns and values 0 to \code{v-1} with coverage
#' strength \code{t} and the number of runs corresponding to the construction
#' (see Details section)
#'
#' @section Details: The function yields designs with columns in v levels based on
#' prime or prime power q as follows:
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
#' @examples
#' ## construction 1
#' ## 2-level
#' CA19.3.19.2 <- cyc(19, 3, 2, 19)
#' coverage(CA19.3.19.2, 3)
#' coverage(CA19.3.19.2, 4)
#'
#' # 3-level
#' CA62.3.32.3 <- cyc(32, 3, 3, 31, type="3a")
#' coverage(CA62.3.32.3, 3)
#' coverage(CA62.3.32.3, 4)
#'
################################################################

#' @export
cyc <- function(k, t, v, q, type=NULL, primitive=NULL){
  ## k is the number of columns
  ## t is the requested coverage
  ## t < k
  ## v is the common number of levels
  ## q is the prime or prime power on which the construction is based
  ## type is the construction type ("1","2","3","3a","3b","4a","4b")

  ## primedat is in sysdata.rda
  qs <- primedat$q
  firstprimitive <- primedat$primitive
  ## define minus for a gf
  gf_minus <- function(x,y,gf){
    ## calculates x - y, after reducing both by a mod operation to
    ## gf entries
    q <- gf$q; p <- gf$p
    x <- x%%q; y <- y%%q
    yn <- gf$neg[y+1]
    gf$plus[x+1,yn+1]
  }
  gf <- lhs::create_galois_field(q)
  ## start vector is the vector of "logarithms of 1:q" (mod v) w.r.t. the
  ##     base chosen as a primitive element (omega) of the group based on q
  if (is.null(type)) type <- "1"
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
  for (i in 1:q)
    for (j in 1:q)
      A[i,j] <- xstart[gf_minus(j,i, gf)+1]
      ## the complicated way of calculating j-i is needed for prime powers,
      ## otherwise (j-i)%%q would do
      ## might be worthwhile to separate, because the %% variant is faster
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
                        function(obj) (rbind(A,0)+obj)%%v)))
  if (type=="4a") return(do.call(rbind, lapply(0:(v-1),
                        function(obj) (cbind(rbind(A,0),0)+obj)%%v)))
  if (type=="4b") return(do.call(rbind, lapply(0:(v-1),
                        function(obj) (cbind(rbind(A,
                                                   matrix(1:(v-1),nrow=v-1,ncol=q)),
                                             1)+obj)%%v)))
}
