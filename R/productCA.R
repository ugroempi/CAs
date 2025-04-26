#' Product constructions
#'
#' Functions to construct a CA from a product construction of two ingoing CAs
#'
#' @rdname productCA
#'
#' @aliases productCA1
#' @aliases productCA2
#' @aliases productPCA
#' @aliases SCA_Bose
#' @aliases CA_to_PCA
#'
#' @usage productCA1(D1, D2, check=TRUE, checklevels=TRUE, dupremove=TRUE, ...)
#' @usage productCA2(A, B, t, check=FALSE, ...)
#' @usage productPCA(D1, D2=NULL, k1=NULL, l1=NULL, start0=TRUE, ...)
#' @usage SCA_Bose(q, ...)
#' @usage CA_to_PCA(D, ...)
#'
#' @param D1 an N x k CA of strength 2 with v levels, split into k1 + k2 columns for \code{productPCA}
#' @param D2 an N x l CA of strength 2 with v levels, split into l1 + l2 columns for \code{productPCA}
#' @param D an N x k CA of strength 2 with v levels
#' @param A  an N x k CA of strength \code{t} with v levels,
#' @param B  an M x k CA of strength \code{t} with w levels,
#' @param t strength of \code{A} and \code{B}
#' @param check logical, if TRUE, checks required strength of ingoing CAs (may substantially increase run time for larger CAs, especially for \code{t}>2)
#' @param checklevels logical, if TRUE, checks that levels of \code{D1} and \code{D2} coincide; this check may be switched off, in order to support interim steps for other constructions
#' @param dupremove logical, if TRUE, removes duplicated rows
#' @param k1 width of P block in \code{D1} (see Details section)
#' @param l1 width of P block in \code{D2} (see Details section)
#' @param start0 logical: do the integer levels start with 0 or with 1?
#' @param q a prime power for the Bose construction
#' @param ... currently not used
#'
#' @section Details:
#' \bold{I don't like the function names, they will likely change.}
#'
#' Function \code{productCA1} yields a CA(N+M, 2, kl, v) from a product construction, as, e.g., described in Torres-Jimènez et al. (2019).\cr
#' Function \code{productCA2} yields a CA(N*M, t, k, v*w) by combining the corresponding columns of the two CAs into \code{v*w} levels, according to Theorem 3.2 of Chateauneuf, Colbourn and Kreher (.\cr
#'
#' @returns Function \code{productCA1} returns a CA(N+M, 2, k*l, v) (N+M x k*l matrix).\cr
#' Function \code{productCA2} yields a CA(N*M, t, k, v*w) (N*M x k matrix).\cr
#' The returned matrices have their levels coded as the ingoing matrices, i.e., integers starting with 0 or 1.
#'
#' Function \code{CA_to_PCA} returns an equivalent CA with columns, rows and levels swapped such that\cr
#' the matrix is partitioned into a PCA structure (Colbourn et al. 2006), arranged here
#' (different from Colbourn et al.) such that\cr
#' the \code{v x k1} matrix P (with each column consisting of a permutation of all levels) makes up the top left block,\cr
#' the top right \code{v x k2} block is constant with all elements equal to the smallest level (for an SCA, which is a special case of a PCA) or has arbitrary levels (for a PCA),\cr
#' and the remaining blocks may have arbitrary elements.\cr
#' Function\code{SCA_Bose} rearranges the columns of a Bose array from \code{\link[lhs]{createBose}} into the column order needed for this package and adds the \code{type} attribute (see below).
#'
#' Function {\code{productPCA}} returns a new PCA in \code{k1*l1 + k1*l2 + k2*l1} columns and \code{N + M - 2} rows.\cr
#'
#' The results of functions \code{CA_to_PCA}, \code{SCA_Bose} and {\code{prodPCA}} have an attribute \code{type}, which is a list with elements \code{type} (PCA or SCA), \code{k1} and \code{k2}.
#'
#'
#' @references Chateauneuf and Kreher (2002), Colbourn et al. (2006), Torres-Jimenez et al. (2019)
#'
#' @examples
#' # product CA (N+M x k*l matrix)
#' A <- cyc(19,2)
#' B <- KSK(k=12)
#' dim(A); dim(B)
#' E <- productCA1(A, B)
#' dim(E)
#' coverage(E, 2)
#' ## more than twice as large as optimal
#' eCAN(2,228,2)
#'
#' ## two CAs with four columns each and strength 3
#' A <- lhs::createBush(3, 4, bRandom=FALSE)  ## strength 3, 4 columns in 3 levels
#' dim(A)
#' B <- lhs::createBush(4, 4, bRandom=FALSE)   ## strength 3, 4 columns in 4 levels
#' dim(B)
#' ## product array with 3*4=12 levels and 27*64=1728 runs
#' E <- productCA2(A, B, 3)
#' dim(E)
#' tail(E)
#' coverage(E,3)
#' eCAN(3,4,12)   ## the CA is optimal (as of April 2025)
#'
#' ## utility CA_to_PCA
#' SCA61.2..9.1..7 <- CA_to_PCA(CS_MS(10, 7))
#' attributes(SCA61.2..9.1..7)
#'
#' CA_to_PCA(DoE.base::L36.2.8.6.3[,-(1:8)])
#'
#' ## utility SCA_Bose
#' SCA_Bose(3)
#'
#' ## apply PCA construction
#' productPCA(SCA_Bose(3))
#'
#' ## for non-prime power case
#' ## use optimum CA(46,2,9,6)
#' MS6 <- CS_MS(k=9, v=6)
#' MS6 <- CA_to_PCA(MS6)
#' head(MS6)  ## it is an SCA with k1=8 and k2=9
#' dim(CA.86.6.80 <- productPCA(MS6, MS6))
#' head(DoE.base::L36.2.8.6.3[,-(1:8)]) ## it is an SCA with k1=3 and k2=0
#' eCAN(2,80,6)
#' dim(productPCA(MS6, DoE.base::L36.2.8.6.3[,-(1:8)] - 1, l1=3))
#' eCAN(2,27,6)
#'
#' @export
productCA1 <- function(D1, D2, check=TRUE, checklevels=TRUE, dupremove=TRUE, ...){
   stopifnot(is.matrix(D1)); stopifnot(is.matrix(D2))
   levs1 <- unique(levels.no(D1))
   levs2 <- unique(levels.no(D2))
   stopifnot(length(levs1)==1);stopifnot(length(levs2)==1)
   if (checklevels){
     stopifnot(levs1==levs2)
     stopifnot(all(unique(c(D1))==unique(c(D2))))
   }
   if(check){
     stopifnot(all(coverage(D1, 2)==1))
     stopifnot(all(coverage(D2, 2)==1))
   }
   k <- ncol(D1); l <- ncol(D2)
   aus <- rbind(D1[,rep(1:k, times=l)],
                D2[,rep(1:l, each=k)])
   if (dupremove) return(aus[!duplicated(aus),])
   aus
}

#' @export
productCA2 <- function(A, B, t, check=FALSE, ...){
   stopifnot(is.matrix(A)); stopifnot(is.matrix(B))
   stopifnot(is.numeric(A)); stopifnot(is.numeric(B))
   stopifnot(min(A)>=0); stopifnot(min(B)>=0)
   stopifnot(min(A)%%1==0); stopifnot(min(B)%%1==0)
   start0 <- TRUE
   if (min(A)==1){
     if(!min(B)==1) stop("Levels of A and B must both start with 0 or both start with 1")
     start0 <- FALSE
   } else
   if(!min(B)==0) stop("Levels of A and B must both start with 0 or both start with 1")
   v <- unique(levels.no(A))
   w <- unique(levels.no(B))
   stopifnot(length(v)==1); stopifnot(length(w)==1)
   k <- ncol(A); l <- ncol(B); stopifnot(k==l)
   N <- nrow(A); M <- nrow(B)
   if(check){
     stopifnot(all(coverage(A, t)==1))
     stopifnot(all(coverage(B, t)==1))
   }
   aus <- matrix(NA, N*M, k)
   if (!start0) {
     A <- A - 1
     B <- B - 1
   }
   for (i in 1:k)
      aus[,i] <- rep(A[,i], each=M)*w + rep(B[,i], v)
   if (!start0) aus <- aus+1
   aus
}

#' @export
CA_to_PCA <- function(D, ...){
  stopifnot(is.matrix(D))
  stopifnot(is.numeric(D))
  v <- levels.no(D)[1]
  k <- ncol(D)
  N <- nrow(D)
  hilf <- maxconstant(D, verbose=2)

  ## entire matrix has v constant rows, i.e., k2=0
  if (length(unique(hilf[v,]))==1){
     aus <- hilf
     aus[1:v,] <- aus[ord(aus[1:v,]),]
     attr(aus, "type") <- list("SCA", k1=k, k2=0)
     return(aus)
  }

  ## not all columns can be made constant for v rows
  ## more constant rows?
    xcol <- NULL

    ## k2 = 1
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
       attr(aus, "type") <- list(type=type, k1=k-1, k2=1)
       return(aus)
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
           attr(aus, "type") <- list(type=type, k1=k-2, k2=2)
           return(aus)
         }else{
           message("D is not a PCA with k2 up to 2. Larger k2 are not yet implemented")
         }
    }
}

#' @export
SCA_Bose <- function(q, ...){
  aus <- lhs::createBose(q, q+1, bRandom=FALSE)[,c(2:(q+1),1)]
  attr(aus, "type") <- list(type="SCA", k1=q, k2=1, eta=3)
  aus
}

#' @export
productPCA <- function(D1, D2=NULL, k1=NULL, l1=NULL, start0=TRUE, ...){
  hilf <- matcheck(D1, start0=start0, PCAcheck = TRUE)
  v <- hilf$v
  k <- hilf$k
  N <- hilf$N
  if (!(identical(hilf$PCAstatus,"unavailable") || is.null(hilf$PCAstatus))){
    type1 <- hilf$PCAstatus$type
    k1_stored <- hilf$PCAstatus$k1
    if (!is.null(k1) && !k1==k1_stored) stop("conflicting information on k1") else
      k1 <- k1_stored
  }else{
    if (is.null(k1)) stop("As D1 does not contain information on its PCA structure,\n k1 must be provided")
    type1 <- "PCA"
    if (k==k1 || all(D1[1:v,(k1+1):k]==min(D1))) type1 <- "SCA"
  }
  if (!is.null(D2)){
    hilf <- matcheck(D2, start0=start0, PCAcheck = TRUE)
    stopifnot(hilf$v==v)
    l <- hilf$k
    M <- hilf$N
    if (!(identical(hilf$PCAstatus,"unavailable") || is.null(hilf$PCAstatus))){
      type2 <- hilf$PCAstatus$type
      if (!type2=="SCA") stop("D2 must be an SCA")
      l1_stored <- hilf$PCAstatus$k1
      if (!is.null(l1) && !l1==l1_stored) stop("conflicting information on l1") else
        l1 <- l1_stored
    }else{
      if (is.null(l1)) stop("As D2 does not contain information on its SCA structure,\n l1 must be provided")
    }
  }else{
    ## no D2 specified, use lhs::createBose
    D2 <- try(SCA_Bose(v))  ## throws error for non-prime v
    if ("try-error" %in% class(D2)) stop("For non-prime number of levels v, D2 must be specified.")
    l <- v+1; M <- v^2; l1 <- v
  }
  ### we now have v, k, l, N, M, k1, l1
  A1 <- D1[(v+1):N,1:k1, drop=FALSE];
  if (k1<k){
      A2 <- D1[(v+1):N,(k1+1):k, drop=FALSE]
      X <- D1[1:v,(k1+1):k, drop=FALSE]
    }else{
      ## empty matrices in the two right-hand side blocks
      A2 <- matrix(1, N-v, 0)
      X <- matrix(1, N-v, 0)
    }
  B1 <- D2[(v+1):M,1:l1, drop=FALSE];
  if (l1 < l){
    B2 <- D2[(v+1):M,(l1+1):l, drop=FALSE]
    } else{
      B2 <- matrix(ifelse(start0,0,1), M-v, 0)
    }
  ## first part, assuming k1>0 and l1>0
  block1 <- rbind(
      matrix(1:v-as.numeric(start0), v, k1*l1), ## P
      productCA1(A1,B1, check=FALSE, checklevels = FALSE, dupremove=FALSE)
    )
  if (k==k1) block2 <- matrix(0, N+M-v, 0) else
    block2 <- rbind(
      do.call("cbind", rep(list(X), l1)),
      productCA1(A2,B1, check=FALSE, checklevels = FALSE, dupremove=FALSE)
    )
  if (l==l1) block3 <- matrix(0, N+M-v, 0) else
    block3 <- rbind(
      matrix(ifelse(start0,0,1), v, (l-l1)*k1),
      productCA1(A1,B2, check=FALSE, checklevels = FALSE, dupremove=FALSE)
    )

  aus <- cbind(block1, block2, block3)
  attr(aus, "type") <- list(type=type1, k1=k1*l1, k2=k1*(l-l1)+(k-k1)*l1)
  aus
}
