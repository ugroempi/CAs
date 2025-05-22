#' Functions for partitioned CAs and their products
#'
#' Functions for partitioned CAs (PCAs) of coverage strength 2, and for a product of two PCAs.
#' The number of rows of a product of a PCA(N, 2, k, v) and a PCA(M, 2, l, v) is N + M - v,
#' the number of columns depends on the splits k = k1+k2 and l = l1+l2.
#'
#' @rdname productPCA
#'
#' @aliases productPCA
#' @aliases SCA_Bose
#' @aliases CA_to_PCA
#' @aliases is.PCA
#'
#' @usage productPCA(D1, D2=NULL, k1=NULL, l1=NULL, ...)
#' @usage SCA_Bose(q, ...)
#' @usage CA_to_PCA(D, tryhard=FALSE, ...)
#' @usage is.PCA(D, start0=NULL, ...)
#'
#' @param D1 an N x k CA of strength 2 with v levels, split into k1 + k2 columns for \code{productPCA}
#' @param D2 an M x l CA of strength 2 with v levels, split into l1 + l2 columns for \code{productPCA}
#' @param D an N x k CA of strength 2 with v levels
#' @param k1 width of P block in \code{D1} (see Details section)
#' @param l1 width of P block in \code{D2} (see Details section)
#' @param q a prime power for the Bose construction
#' @param tryhard logical: if \code{TRUE}, \code{CA_to_PCA} tries to achieve \code{v} constant rows for \code{k-1} or \code{k-2} columns; this may take a while in case of many columns.
#' @param start0 \code{NULL}, or logical: do the values start with 0 (otherwise with 1)?
#' @param ... currently not used
#'
#' @section Details:
#' Function \code{productPCA} exploits a partitioned structure of the two CAs it combines
#' (see Section "Partitioned Covering Array (PCA)"). The construction can be viewed as a kind of product,
#' but is sometimes also called "cut-and-paste". It was proposed in Colbourn et al. (2006)
#' and was slightly generalized later, see, e.g., in Colbourn and Torres-Jimenez (2013).
#'
#' Function \code{productPCA} also serves as work horse function for combining two CAs in many instances of
#' the \code{\link{CAEX}} construction.
#'
#' \code{CA_to_PCA} rearranges \code{D} into a PCA structure (see Section "Partitioned Covering Array (PCA)")
#' by rearranging columns, swapping symbols, and in case of \code{tryhard=TRUE} also rearranging rows.\cr
#' In case of \code{tryhard=TRUE}, it is tried to achieve PCA status \code{k2=1} or \code{k2=2}; in that case,
#' if it is possible to achieve \code{v} constant rows (i.e., \code{k2=0}), the function finds these.\cr
#' Without \code{tryhard=TRUE}, the resulting PCA will have \code{k1} as the number of columns with \code{v} different
#' levels in the first \code{v} rows.\cr
#' \code{tryhard=TRUE} may be useful for small designs, but may substantially increase run times for \code{D} with many rows or many columns.
#'
#' Function \code{SCA_Bose} obtains the strength 2 Bose constructed CA (q+1 columns in q levels, q^2 runs, q a prime power)
#' and rearranges it into SCA shape.
#'
#' Function \code{CA_to_PCA} brings its argument into the best possible PCA shape by swapping columns only.
#' If the argument \code{tryhard} is invoked, the function tries to achieve \code{k1=k-1} or \code{k1=k-2} by trying \code{\link{maxconstant}} on
#' subsets of columns. Function \code{is.PCA} checks the PCA status of its argument without modifying it in any way.
#'
#' @section Partitioned Covering Array (PCA):
#' For an \code{N x k} CA(N, t, k, v) to be a PCA, the number of columns is split into \code{k = k1 + k2}, with \code{k1>0} (otherwise no PCA),
#' the number of rows is split into \code{N = v + (N - v)}; such a PCA is denoted as PCA(N, t, (k1,k2), v).\cr
#' The PCA structure, in this package, is arranged with top and bottom parts swapped versus
#' partitioning of Colbourn et al. (2006). Parts are arranged here such that\cr
#' the \code{v x k1} matrix P (with each column consisting of the identity permutation of all levels)
#' makes up the top left block,\cr
#' the top right \code{v x k2} block (which may be empty) is constant with all elements equal to the smallest level (for an SCA, which is a special case of a PCA)
#' or has arbitrary levels (for a PCA),\cr
#' and the remaining two blocks (bottom left and bottom right) may arbitrary elements.
#'
#' @returns
#' Function {\code{productPCA}} returns a new PCA in \code{k1*l1 + k1*l2 + k2*l1} columns and \code{N + M - v} rows
#' (matrix of class \code{ca}).\cr
#' The returned matrix has its levels coded like the ingoing matrices, i.e., integers starting with 0 or 1.
#'
#' Function \code{CA_to_PCA} returns an equivalent CA of S3 class \code{ca} that has PCA (or SCA) structure (see Section "Partitioned Covering Array (PCA)").
#'
#' Function\code{SCA_Bose} creates a Bose array from \code{\link[lhs]{createBose}} and rearranges its columns
#' into SCA structure (see Section "Partitioned Covering Array (PCA)").
#'
#' The results of functions \code{CA_to_PCA}, \code{SCA_Bose} and {\code{prodPCA}} have an attribute \code{PCAstatus},
#' which is a list with elements \code{type} (PCA or SCA), \code{k1} and \code{k2}.
#'
#' Function \code{is.PCA} returns a logical with the analogous attribute; it can be (mis)used for checking
#' consistency between the codings of designs by specifying \code{start0}.
#'
#' @references Colbourn et al. (2006), Colbourn and Torres-Jimenez (2013), Torres-Jimenez et al. (2019)
#'
#' @examples
#' #################################################
#' ### productPCA
#' #################################################
#' A <- SCA_Bose(4)  ## Bose array, arranged as SCA
#' EPCA <- productPCA(A, A)
#' coverage(EPCA, 2)
#' dim(EPCA)  ## 24 columns only, 28 runs only
#' eCAN(2, 24, 4) ## optimal, this construction!
#' eCAK(2, 28, 4) ## optimal, this construction
#'
#' ########################################
#' ## utilities CA_to_PCA and is.PCA
#' ########################################
#' SCA61.2..9.1..7 <- CA_to_PCA(CS_MS(10, 7), tryhard=TRUE)
#' attributes(SCA61.2..9.1..7)
#'
#' is.PCA(DoE.base::L36.2.8.6.3[,-(1:8)])
#' CA_to_PCA(DoE.base::L36.2.8.6.3[,-(1:8)])
#'
#' ## utility SCA_Bose
#' SCA_Bose(3)
#'
#' ## apply PCA construction
#' E <- productPCA(SCA_Bose(3))
#' dim(E)
#' coverage(E,2)
#' eCAN(2,15,3)  ## optimal for 15 columns
#' eCAK(2,15,3)  ## but does not cover 20 columns
#'
#' ## for non-prime power case
#' ## use optimum CA(46,2,9,6)
#' MS6 <- CS_MS(k=9, v=6)
#' MS6 <- CA_to_PCA(MS6)
#' MS6 <- CA_to_PCA(MS6, tryhard=TRUE)
#' head(MS6)  ## it is an SCA with k1=8 and k2=9
#' dim(CA.86.6.80 <- productPCA(MS6, MS6))
#' ## four too many runs
#' eCAN(2,80,6)
#'
#' ## another case
#' head(DoE.base::L36.2.8.6.3[,-(1:8)]) ## it is an SCA with k1=3 and k2=0,
#'                                      ## as all columns have six distinct
#'                                      ## values in first six rows
#' dim(productPCA(MS6, CA_to_PCA(DoE.base::L36.2.8.6.3[,-(1:8)] - 1)))
#' eCAN(2,27,6)  ## 8 too many runs
#'
#' ## the function also works for designs with flexible values (which are denoted by NA)
#' ## productPCA
#' D <- productPCA(CAEX(N=16), CAEX(N=16))
#' dim(D)
#' eCAN(2,405,3) ## slightly worse than the best CA, which is obtained as CAEX(N=27)
#' is.PCA(D) ## is already good and cannot be further improved by the simple method used in CA_to_PCA
#' \dontrun{coverage(D, 2)} ## slightly longer run time
#'

#' @export
productPCA <- function(D1, D2=NULL, k1=NULL, l1=NULL, ...){
  Call <- sys.call()
  hilf <- matcheck(D1, PCAcheck = TRUE)
  v <- hilf$v
  k <- hilf$k
  N <- hilf$N
  flexible <- hilf$flexval
  if (!is.null(flexible)){
    if (!is.na(flexible)) D1[D1==flexible] <- NA
    flexible <- NA
  }
  start0 <- hilf$start0
  k1_stored <- try(hilf$PCAstatus$k1) ## read from hilf
  if ("try-error" %in% class(k1_stored))
    stop("D1 does not seem to be a PCA")
  type1 <- hilf$PCAstatus$type
  if (!is.null(k1)){
    if (!k1==k1_stored) stop("conflicting information on k1")
  }else
    k1 <- k1_stored

  ## should not happen any more
  if (is.null(k1)) stop("As D1 does not contain information on its PCA structure,\n k1 must be provided")
  if (is.null(type1)){
    type1 <- "PCA"
    if (k==k1 || all(D1[1:v,(k1+1):k]==min(D1, na.rm=TRUE)))
      type1 <- "SCA"
  }
  ## now handle D2
  if (!is.null(D2)){
    hilf <- matcheck(D2, start0=start0, PCAcheck = TRUE)
    stopifnot(hilf$v==v)
    l <- hilf$k
    M <- hilf$N
    flexval2 <- hilf$flexval
    if (!is.null(flexval2))
      if (!is.na(flexval2)){
        D2[D2==flexval2] <- NA
        flexible <- NA
      }
    l1_stored <- try(hilf$PCAstatus$k1) ## read from hilf
    type2 <- hilf$PCAstatus$type
    if ("try-error" %in% class(l1_stored))
      stop("D2 does not seem to be a PCA")
    type2 <- hilf$PCAstatus$type
    if (!is.null(l1)){
      if (!l1==l1_stored) stop("conflicting information on l1")
    }else
      l1 <- l1_stored

    ## should not happen any more
    if (is.null(l1))
      stop("As D2 does not contain information on its PCA structure,\n l1 must be provided")
  }else{
    ## no D2 specified, use Bose construction (works for prime power v only)
    D2 <- try(SCA_Bose(v))  ## throws error for non-prime power v
    if ("try-error" %in% class(D2))
      stop("For non-prime number of levels v, D2 must be specified.")
    l <- v+1; M <- v^2; l1 <- v; type2 <- "SCA"
  }
  ### we now have v, k, l, N, M, k1, l1, start0, flexible
  ### The construction

  ## left blocks of D1
  P1 <- matrix(1:v-as.numeric(start0), v, k1)
  A1 <- D1[(v+1):N, 1:k1, drop=FALSE];
  if (k1 < k){
    ## lower right block of D1
    A2 <- D1[(v+1):N, (k1+1):k, drop=FALSE]
    ## upper right block of D1
    X <- D1[1:v,(k1+1):k, drop=FALSE]
  }else{
    ## empty matrices in the two right-hand side blocks
    A2 <- matrix(1, N-v, 0)
    X <- matrix(1, v, 0)
  }
  ## left blocks of D2
  P2 <- matrix(1:v-as.numeric(start0), v, l1)
  B1 <- D2[(v+1):M,1:l1, drop=FALSE];
  if (l1 < l){
    ## lower right block of D2
    B2 <- D2[(v+1):M,(l1+1):l, drop=FALSE]
    ## upper right block of D2
    Y <- D2[1:v,(l1+1):l, drop=FALSE]
  } else{
    ## empty matrices on the two right-hand side blocks
    B2 <- matrix(1, M-v, 0)
    Y <- matrix(1, v, 0)
  }
  ## first part, assuming k1>0 and l1>0
  ## left block, width k1*l1
  block1 <- rbind(
    matrix(1:v-as.numeric(start0), v, k1*l1), ## P
    productCA_raw(A1, B1)
  )
  ## middle block
  if (k==k1) block2 <- matrix(0, N+M-v, 0) else
    block2 <- rbind(
      do.call("cbind", rep(list(X), l1)),
      productCA_raw(A2, B1)
    )
  ## right block
  if (l==l1) block3 <- matrix(0, N+M-v, 0) else
    block3 <- rbind(
      kronecker(Y, matrix(1,1,k1)),
      productCA_raw(A1, B2)
    )

  aus <- cbind(block1, block2, block3)
  attr(aus, "PCAstatus") <- list(type=type1, k1=k1*l1, k2=k1*(l-l1)+(k-k1)*l1)
  class(aus) <- c("ca", class(aus))
  attr(aus, "Call") <- Call
  aus
}

#' @export
is.PCA <- function(D, start0=NULL, ...){
  hilf <- matcheck(D, start0=start0)
  v <- hilf$v
  k <- hilf$k
  N <- hilf$N
  start0 <- hilf$start0
  flexible <- hilf$flexvalue
  # ### remove, because this does not make sense in this function
  # if (!is.null(flexible)){
  #   if (is.na(flexible)) profile <- apply(D, 2, function(obj) sum(is.na(obj))) else
  #   profile <- apply(D, 2, function(obj) sum(obj == flexible))
  # }
  ## check whether the file is already in PCA / SCA form
  comparevec <- 1:v - as.numeric(start0)

  ## increase k1 for each column found, from left to right
  k1 <- 0
  for (i in 1:k){
    if (all(D[1:v, i]==comparevec)) k1 <- k1+1 else break
  }
  k2 <- k - k1
  if (k1 > 1){
    type <- "PCA"
    if (k1 < k){
      if (all(D[1:v, (k1+1):k]==1-as.numeric(start0))) type <- "SCA"
    }
    aus <- TRUE
  }else aus <- FALSE
  ## return TRUE with info attached
  if (aus) attr(aus, "PCAstatus") <- list(type=type, k1=k1, k2=k2) else attr(aus, "PCAstatus") <- FALSE
  # if (!is.null(flexible)) attr(aus, "flexible") <- list(value=flexible, profile=profile)
  aus
}

#' @export
CA_to_PCA <- function(D, tryhard=FALSE, ...){
  Call <- sys.call()
  ## assumes that undeclared flexible values are coded as NA
  hilf <- matcheck(D, flexible=NA)
  v <- hilf$v; k <- hilf$k; N <- hilf$N; start0 <- hilf$start0; flexible <- hilf$flexval
  hilf <- maxconstant(D, verbose=2)

  ## entire matrix has v constant rows, i.e., k2=0
  if (length(unique(hilf[v,]))==1){
    aus <- hilf
    aus[1:v,] <- aus[ord(aus[1:v,]),]
    attr(aus, "PCAstatus") <- list("SCA", k1=k, k2=0)
    return(aus)
  }

  ## determine the column indices for which the first row has
  ## v distinct values
  hilf <- levels.no(D[1:v,])
  good <- which(hilf==v)
  kcheap <- length(good)

  ## there is already a PCAstatus attribute
  if (!is.null(PCAstatus <- attr(D, "PCAstatus"))){
    if (!is.list(PCAstatus))
      attr(D, "PCAstatus") <- "unknown" else{
      if (PCAstatus$type %in% c("PCA","SCA")){
        ## valid attribute is assumed to be correct
        if (PCAstatus$k1>=kcheap){
          message("The properties of the PCA D cannot be improved by function CA_to_PCA.")
          return(D)  ## return unchanged input matrix
        }
      }else attr(D, "PCAstatus") <- "unknown"
      }
    }

  ## with valid non-NULL attribute PCAstatus that is as good as
  ##    possible without row swaps, the array was already returned
  ## with any other non-NULL attribute PCAstatus or FALSE,
  ##    the attribute was set to "unknown"

    ## run check after moving the columns
    ## with distinct values in the first v rows
    ## to the left
    D <- D[,c(good, setdiff(1:k, good))]
    if (kcheap > 0){
      for (i in 1:kcheap){
          D <- permvals(D, i, D[1:v, i], 1:v - as.numeric(start0))
      }
    PCAstatus <- is.PCA(D, start0=start0, flexible=NA)
    attr(D, "PCAstatus") <- attr(PCAstatus, "PCAstatus")
    attr(D, "flexible") <- attr(PCAstatus, "flexible")
    }else{
      if (!tryhard) {
        message("The first ", v, " rows do not have any columns that can serve as the k1-part of a PCA.")
        attr(D, "PCAstatus") <- FALSE
      }
    }
    if (!tryhard) return(D)


  ## not all columns can be made constant for v rows
  ## and there is no valid PCAstatus attribute

  ## attempt to make k-1 or k-2 columns constant in v rows
  ## more constant rows?

  ## tryhard=TRUE
  xcol <- NULL

  ## try to achieve a split with k1=k-1 and k2=1
  for (i in 1:k){
    hilf <-  maxconstant(D[,-i], verbose=2)
    if (length(unique(hilf[v,]))==1){
      xcol <- i
      break
    }
  }

  if (!is.null(xcol)){
    aus <- D[,c(setdiff(1:k,xcol),xcol)]
    ## constantrow-indices used for hilf
    inds <- attr(hilf, "constant_rows")$row_set
    aus <- rbind(cbind(hilf[1:v,], D[inds,xcol,drop=FALSE]),
                 cbind(hilf[-(1:v),], D[setdiff(1:N, inds), xcol, drop=FALSE]))
    k1 <- k-1; k2 <- 1
    aus[1:v,] <- aus[ord(aus[1:v,]),]  ## sort into entries 1:v for first k1 columns
    type <- "PCA"
    if (length(unique(aus[1:v,(k1+1):(k1+k2)]))==1){
      type <- "SCA"
      if (!aus[1,k1+1]==min(aus)){
        from <- aus[1,k1+1]; to <- min(aus)
        for (i in 1:k2){
          aus <- swapvals(aus, k1+i, from, to)
        }
      }
    }
    attr(aus, "PCAstatus") <- list(type=type, k1=k-1, k2=1)
  } else{
    ## try pairs of columns
    pickpairs <- nchoosek(k,2)
    for (i in 1:ncol(pickpairs)){
      hilf <-  maxconstant(D[,-pickpairs[,i]], verbose=2)
      if (length(unique(hilf[v,]))==1){
        xcol <- pickpairs[,i]
        break ## breaks the loop over i
      }
    }
    if (!is.null(xcol)){
      aus <- D[,c(setdiff(1:k,xcol),xcol)]
      ## constantrow-indices used for hilf
      inds <- attr(hilf, "constant_rows")$row_set
      aus <- rbind(cbind(hilf[1:v,], D[inds,xcol,drop=FALSE]),
                   cbind(hilf[-(1:v),], D[setdiff(1:N, inds), xcol, drop=FALSE]))

      k1 <- k-1; k2 <- 1
      aus[1:v,] <- aus[ord(aus[1:v,]),]  ## sort into entries 1:v for first k1 columns
      type <- "PCA"
      if (length(unique(aus[1:v,(k1+1):(k1+k2)]))==1){
        type <- "SCA"
        if (!aus[1,k1+1]==min(aus)){
          from <- aus[1,k1+1]; to <- min(aus)
          for (i in 1:k2){
            aus <- swapvals(aus, k1+i, from, to)
          }
        }
      }
      attr(aus, "PCAstatus") <- list(type=type, k1=k-2, k2=2)
    }else{
      message("D is not a PCA with k2 up to 2. tryhard did not improve the PCA status.")
    }
  }
  attr(aus, "Call") <- c(attr(aus, "Call"), Call)
   class(aus) <- c("ca", class(aus))
   aus
}

#' @export
SCA_Bose <- function(q, ...){
  Call <- sys.call()
  aus <- lhs::createBose(q, q+1, bRandom=FALSE)[,c(2:(q+1),1)]
  attr(aus, "PCAstatus") <- list(type="SCA", k1=q, k2=1)
  attr(aus, "Call") <- Call
  class(aus) <- c("ca", class(aus))
  aus
}
