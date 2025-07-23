#' Function to create a strength 4, 5, or 6 CA for 2-level columns based on Colbourn-Keri cover starters
#'
#' A strength 4, 5, or 6 CA in k 2-level columns is created
#' from a starter vector of length k or k-1. It has 2*k or 2*(k-1) runs.
#' The construction uses cyclic permutations and is based on Colbourn and Keri (2009).
#'
#' @rdname CS_CK
#'
#' @aliases CS_CK
#' @aliases N_CS_CK
#' @aliases k_CS_CK
#'
#' @usage CS_CK(k, t=4, v=2, starter=NULL, method=NULL, start0=TRUE, ...)
#' @usage N_CS_CK(t, k, v=2)
#' @usage k_CS_CK(t, N, v=2)
#'
#' @param k the number of columns of the output CA
#' @param t strength of the output CA, can be 4, 5, or 6
#' @param v the number of levels of the output CA (always 2)
#' @param starter \code{NULL} for automatic choice (recommended), or a vector of length at least \code{k-1};\cr
#'     if a vector is provided, it is the users responsibility to ensure the desired properties of
#'     the result; the specification of a value for \code{t} has no effect whatsoever in that case
#' @param method \code{NULL}, or the specification of how to process a custom starter vector;
#'     if non-NULL, can be either a number from 1, 2, 3 or a character string from
#'     "CS_CK1", "CS_CK2", "CS_CK3";\cr
#'     \code{method} is ignored, unless \code{starter} is specified, in which case it is required
#' @param start0 logical; if TRUE, values of the output CA start at 0, otherwise at 1;
#'     values of the starter must always start at 0
#' @param N affordable run size
#' @param ... currently not used
#'
#' @details
#' The function implements the designs provided in Colbourn and Keri (2009).
#' The \code{method} numbers refer to
#' the lemmata in that paper: Methods 1 and 3 yield \code{N=2*k},
#' method 2 yields \code{N=2*(k-1)};
#' the starter vector has length \code{k} for method 1
#' and length \code{k-1} for methods 2 and 3.
#' Method 3 always yields designs with the last two rows constant.
#'
#' Besides methods 1, 2 and 3, the function also provides constructions for \code{t=4,5} by derivation
#' from strengths 5 and 6, respectively (method \code{derive} in the data.frame object \code{ColbournKeriCombis},
#' that provides an overview of the constructions).
#' The cover starters for strength 4 are in the object \code{ColbournKeriStarters},
#' the cover starters for higher strengths are calculated with the internal function \code{cycvec} (that is
#' also used for the cyclotomy constructions of function \code{\link{cyc}}).
#'
#' @returns \code{CS_CK} returns a matrix of class \code{ca} with attributes.
#'
#' @references Colbourn and Keri (2009)
#'
#' @examples
#' eCAN(4,21,2)  ## this is the best construction, it uses method 1
#' Ns(4,21,2)
#' D <- CS_CK(21)  ## strength 4 is the default
#' dim(D)
#' # coverage(D, 4)
#'
#' Ns(4, 40, 2)
#' dim(CS_CK(40)) ## almost as good as best
#' ## construction was derived, stems from a strength 5 CA made with method 3
#'
#' Ns(6, 430, 2)
#' dim(CS_CK(430, t=6)) ## best, uses method 2
#' ## closely related to the Paley type construction by Colbourn
#'
#' ## a few Colbourn Keri-related entries in the Colbourn tables apparently lack
#' ## a postop suffix (status: June 2025)
#' eCAN(5, 359, 2)
#' eCAN(6, 360, 2)
#' ## deriving the 720 run design yields 360 runs
#' ## the 359 runs of the Colbourn tables must be the result of a postop step
#' Ns(5, 35, 2)
#' Ns(5, 359, 2)
#' Ns(5, 503, 2)  ##  from deriving 504 columns in 1008 runs
#'

CS_CK1 <- function(k, starter, ...){
  ## internal function, no input checks needed
  ## works for k<=26, with constructions for k=21, 24, 25, 26
  ## starter vector has length k
  ## Lemma 1 of Colbourn and Keri (2009)
  ## yields strength 4
  aus <- circular(starter)
  rbind(aus, 1-aus)
}

CS_CK2 <- function(k, starter, ...){
  ## internal function, no input checks needed
  ## works for k<=34, with constructions for k=24, 27 to 34
  ## starter vector has length k-1
  ## Lemma 2 of Colbourn and Keri (2009)
  ## yields strength 4

  aus <- cbind(circular(starter))
  cbind(rbind(aus, 1-aus),rep(0:1, each=k-1))
}

CS_CK3 <- function(k, starter, ...){
  ## internal function, no input checks needed
  ## works for k<=27, with constructions for k=12, 20, 22, 23, 25 to 27
  ## starter vector has length k-1
  ## Lemma 3 of Colbourn and Keri (2009)
  aus <- CS_CK2(k, starter=starter)
  rbind(aus, 0, 1)
}

#' @export
CS_CK <- function(k, t=4, v=2, starter=NULL, method=NULL, start0=TRUE, ...){
  Call <- sys.call()
  stopifnot(is.numeric(t), is.numeric(k), is.numeric(v))
  if (!v==2) stop("CS_CK is for v=2 only")
  if (!(k%%1==0)) stop("k must be an integer number")
  stopifnot(t %in% 4:6)
  if (!is.null(starter) && is.null(method)) stop("if starter is given, method must also be given")
  if (is.null(starter)){
    if (!is.null(method)) message("if starter is NULL, method is ignored")
    method <- NULL
  }
  if (!is.null(starter)){
    stopifnot(all(starter %in% c(0,1)))
    stopifnot(length(starter) >= k-1)
    stopifnot(method %in% 1:3 || method %in% c("CS_CK1", "CS_CK2", "CS_CK3"))
     ## 1: CS_CK1 q=k, N=2*q=2*k
     ## 2: CS_CK2 q=k-1, N=2*q=2*(k-1)
     ## 3: CS_CK3 q=k-1, N=2*(q+1)=2*k
    if (method %in% 1:3) method <- paste0("CS_CK", method)
    message("With a specified starter, the user is responsible for ensuring coverage.\nSpecifying t does not do it for you.")
    aus <- switch(method, NA, CS_CK1=CS_CK1(k, starter=starter),
                              CS_CK2=CS_CK2(k, starter=starter),
                              CS_CK3=CS_CK3(k, starter=starter))
    if (!start0) aus <- aus+1
    class(aus) <- c("ca", class(aus))
    attr(aus, "Call") <- Call
    attr(aus, "origin") <- paste0("Colbourn Keri, custom starter, method ", method)
    attr(aus, "starter") <- starter
    return(aus[,1:k])
  }
  ## now for built-in starter
    zeilen <- ColbournKeriCombis[ColbournKeriCombis$t>=t & ColbournKeriCombis$k>=k,]
    ## lower strength is at the top,
    ## higher strength at the bottom
    zrows <- nrow(zeilen)
    if (zrows == 0) stop("The request cannot be served by function CS_CK.")
    zeile <- zeilen[zrows + 1 - which.min(rev(zeilen$N)), ,drop=FALSE]
    aus <- eval(parse(text = zeile$code))
    if (!start0) aus <- aus + 1
    class(aus) <- c("ca", class(aus))
    attr(aus, "Call") <- Call
    attr(aus, "origin") <- paste0("Colbourn Keri, strength t=", t, ", ")
    attr(aus, "t") <- t
    aus[,1:k]
}

#' @export
N_CS_CK <- function(t,k,v=2){
  hilf <- ColbournKeriCombis[ColbournKeriCombis[,"t"]==t & ColbournKeriCombis[,"k"]>=k &
                               ColbournKeriCombis[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(min(hilf[,"N"]))
}

#' @export
k_CS_CK <- function(t,N,v=2){
  hilf <- ColbournKeriCombis[ColbournKeriCombis[,"t"]==t & ColbournKeriCombis[,"N"]<=N &
                               ColbournKeriCombis[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  else return(max(hilf[,"k"]))
}


