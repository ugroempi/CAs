#' function to create a strength 2 CA based on recursive multiplication of Bose CAs
#'
#' recBoseCA applies productCA or productPCA recursively to Bose CAs. The number of levels must be a prime or prime power.
#'
#' @rdname recursiveBose
#'
#' @aliases recBoseCA
#' @aliases recursiveBose
#' @aliases N_k_recursiveBose
#' @aliases N_d_recursiveBose
#'
#' @usage recBoseCA(t=2, k=NULL, v=4, type="PCA", ...)
#' @usage recursiveBose(q, d=NULL, k=NULL, type="PCA", ...)
#' @usage N_k_recursiveBose(q, d=2, type="PCA", ...)
#' @usage N_d_recursiveBose(q, k, type="PCA", ...)
#'
#' @param t the requested strength; must be 2
#' @param k \code{NULL}, or the number of columns;\cr
#'          if \code{NULL}, the maximum possible number for \code{d} is returned
#'          (with \code{d} defaulting to 2, if not specified).
#' @param v the same as \code{q}
#' @param q a prime or prime power number of levels
#' @param type \code{"CA"} or \code{"PCA"},
#'          indicating whether \code{\link{productCA}} or
#'          \code{\link{productPCA}} should be used for the product
#' @param d \code{NULL}, or a small positive integer;\cr
#'          \code{d=1} returns \code{SCA_Bose(q)} (q^2 x (q+1)),
#'          \code{d} > 1 the product (\code{\link{productCA}} or \code{\link{productPCA}}, depending on \code{type})
#'          of \code{d} \code{SCA_Bose(q)};\cr
#'          For \code{d=NULL}, the minimum \code{d} for accommodating \code{k} is used.\cr
#'          If \code{k} is also \code{NULL}, the default is \code{d=2}.
#' @param ... currently not used
#'
#' @section Details:
#' The function \code{recBoseCA} uses function \code{recursiveBose}, which in turn uses
#' \code{\link{SCA_Bose}} for creating the \code{q^2 x (q+1)} ingredient matrix,
#' and either function \code{\link{productPCA}} or \code{\link{productCA}}
#' for obtaining the recursive product (depending on the choice for argument \code{type}).\cr
#' For a given $\code{v}$, the function can construct\cr
#' up to \code{v^d+d*v^(d-1)} columns in \code{d*v^2-(d-1)*v} runs (\code{type=="PCA"})\cr
#' or up to \code{(v+1)^d} columns in up to \code{d*v^2-(d-1)*2} runs (\code{type=="CA"}),\cr
#' where \code{d} is the number of Bose CAs that are combined.
#'
#' @returns \code{recBoseCA} returns a uniform CA strength 2 CA with \code{k} \code{v} level columns
#' (matrix of class \code{ca}); its rowsize can be determined a-priori with function \code{\link{N_recBoseCA}}.\cr
#' \code{recursiveBose} returns a uniform strength 2 CA with \code{q} level columns.\cr
#' For \code{type="PCA"} (default), the result has \code{d*q^2 - (d-1)*q} rows and
#' \code{q^d + d*q^(d-1)} columns, levels are coded from 0 to \code{q - 1}.\cr
#' For \code{type="CA"}, the result has \code{d*q^2 - (d-1)*2} rows (upper bound,
#' may be slightly smaller depending on constant rows) and
#' \code{(q+1)^d} columns, levels are coded from 0 to \code{q - 1}.\cr
#' \code{N_k_recursiveBose} returns a named vector with the number of rows and the number of columns k
#' for given \code{q} and \code{d}.\cr
#' \code{N_d_recursiveBose} returns a named vector with the run size \code{N}, the \code{d} for its creation for the requested number of columns \code{k}, and the maximum possible number of columns \code{kmax} for that construction.\cr
#'
#' @examples
#' D <- recBoseCA(k=12)   ## v=4 is the default, and type="PCA"
#' dim(D)
#' coverage(D,2)
#' eCAN(2, 12, 4)
#' Ns(2, 12, 4) ## the optimum CA is available
#'
#' D <- recBoseCA(k=25)   ## v=4 is the default, and type="PCA"
#' dim(D)
#' coverage(D,2)
#' eCAN(2, 25, 4)
#' Ns(2, 25, 4) ## the optimum CA is available
#' ## type="CA" gets close to optimal here
#' ## because k=25 makes it still work with d=2 SCA_Bose
#' D <- recBoseCA(k=25, type="CA")
#' dim(D)
#'
#' dim(D <- recursiveBose(3))
#' eCAN(2,15,3)  ## optimal
#' dim(recursiveBose(3, d=3))
#' eCAN(2,54,3)  ## reasonably close to optimal
#' dim(recursiveBose(3, d=3, type="CA"))
#' eCAN(2,64,3)  ## reasonably close to optimal
#'
#' ## query functions
#' N_k_recursiveBose(4)      ## d=2
#' N_k_recursiveBose(4, 3)   ## d=3
#' N_d_recursiveBose(4, 12)  ## k=12
#'
#' @export
recBoseCA <- function(t=2, k=NULL, v=4, type="PCA", ...){
  Call <- sys.call()
  stopifnot(t==2)
  if (!(v %in% primedat$q)) stop("v must be prime or prime power")
  stopifnot(type %in% c("CA", "PCA"))
  suppressMessages(
    if (is.null(k))
      k <- unname(N_k_recursiveBose(v, type=type)["k"])
    else
     stopifnot(is.numeric(k), k>0, k%%1==0)
    )
  aus <- recursiveBose(v, k=k, type=type)
  class(aus) <- c("ca", class(aus))
  attr(aus, "Call") <- Call
  attr(aus, "t") <- 2
  aus
}
#'
#' @export
recursiveBose <- function(q, d=NULL, k=NULL, type="PCA", ...){
  Call <- sys.call()
  start0 <- TRUE ## default
  if (!is.null(d)){
    stopifnot(is.numeric(d))
    stopifnot(d>0)
    stopifnot(d %% 1 == 0)
  }
  if (is.null(d) && is.null(k)) d <- 2
  if (is.null(d)){
    ## infer d from k, which is not NULL
    stopifnot(is.numeric(k))
    stopifnot(k>0)
    stopifnot(k %% 1 == 0)
    if (type=="CA") d <- ceiling(log(k, base=q+1)) else{
      hilf <- q^(1:6) + (1:6)*q^((1:6)-1)
      d <- min(which(hilf>=k))
    }
    }
    if (d==1) return(SCA_Bose(q))
    A <- SCA_Bose(q)
    for (i in 2:d){
      if (type=="CA") A <- productCA(A, SCA_Bose(q)) else
        A <- productPCA(A, SCA_Bose(q))
    }
  class(A) <- c("ca", class(A))
  attr(A, "Call") <- Call
  attr(A, "origin") <- c("recursiveBose", type)
  attr(A, "t") <- 2
  A
}
#' @export
N_k_recursiveBose <- function(q, d=2, type="PCA", ...){
  if (!q %in% primedat$q)
    stop("recursiveBose requires q to be a prime or prime power")
  stopifnot(type %in% c("PCA", "CA"))
  ## check for d
  stopifnot(is.numeric(d))
  stopifnot(d > 0)
  stopifnot(d %% 1 == 0)
  if (type=="PCA"){
    return(c(N=d*q^2-(d-1)*q, k=q^d + d*q^(d-1)))
  } else{
    message("N can be slightly smaller because of constant rows")
    return(c(N=d*q^2-(d-1)*2, k=(q+1)^d))
  }
}

#' @export
N_d_recursiveBose <- function(q, k, type="PCA", ...){
  if (!q %in% primedat$q)
    stop("recursiveBose requires q to be a prime or prime power")
  stopifnot(type %in% c("PCA", "CA"))
  ## check for k
  stopifnot(is.numeric(k))
  stopifnot(k>1)
  stopifnot(k %% 1 == 0)
  ## figure out the necessary d
  if (type=="CA") d <- ceiling(log(k, base=q+1)) else{
    d <- 1
    kmax <- q+1
    while(k > kmax){
      d <- d + 1
      kmax <- q^d + d*q^(d-1)
    }
  }
  if (type=="PCA") return(c(N=d*q^2-(d-1)*q, d=d, kmax= q^d + d*q^(d-1))) else{
    message("N can be slightly smaller because of constant rows")
    return(c(N=d*q^2-(d-1)*2, d=d, kmax=(q+1)^d))
  }
}
