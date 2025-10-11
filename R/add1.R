#' function to add a columns
#'
#' increases the number of columns by 1,
#' adding at most bestN(t-2,k-1,v)*v*(v-1) rows.
#'
#' @rdname add1
#'
#' @aliases add1CA
#' @aliases N_add1CA
#' @aliases add1
#'
#' @usage add1CA(t, k, v, ...)
#' @usage N_add1CA(t, k, v)
#' @usage add1(D, t, v, ...)
#'
#' @param t integer: the required strength, or the strength of \code{D} (not checked)
#' @param k integer: the required number of columns
#' @param v integer: the required number of levels, or the number of levels of \code{D}
#' @param D a strength \code{t} CA with \code{v} levels in \code{k-1} columns\cr
#' (integer levels, starting at 0 or 1)
#' @param ... currently not used
#'
#' @returns \code{add1} and \code{add1CA} return a CA of strength \code{t} in \code{v} levels in
#' \code{k} columns.\cr
#' It will have at most \code{\link{bestN}(t-2,k-2,v)*v*(v-1)} additional rows.\cr
#' \code{N_add1CA} returns the number of runs obtainable by \code{add1CA}.
#'
#' @section Details:
#' Function \code{add1CA} applies function \code{add1} the the best available
#' strength \code{t} CA at \code{v} levels with \code{k-1} columns in the role of \code{D}.\cr
#' The construction relies on appending to \code{D} a copy of the last column of \code{D},
#' followed by adding rows that consist of\cr
#' the best strength \code{t-2} CA at \code{v} levels with \code{k-2} columns, repeated \code{v*(v-1)} times, in the first \code{k-2} columns,\cr
#' each combined with one of the \code{v(v-1)} distinct level pairs for the last two columns.\cr
#' For \code{t=2}, the best CA is simply a constant row of flexible values.
#'
#' The construction is not interesting for small CAs, but may be useful, where
#' good large CAs, e.g., ones constructed by cyclotomy, have one less column
#' than needed (and a good CA with strength \code{t-2} for \code{k-1} columns is
#' available).
#'
#' @references CKRS, Theorem 3.2
#'
#' @examples
#' D <- powerCA(3,1331,2)
#' dim(D)
#' D2 <- add1(D, 3, 2)
#' dim(D2)
#' eCAN(3,1331,2)
#' eCAN(3,1332,2)
#' N_add1CA(3,1332,2)
#'
#' dim(add1CA(3,1332,2))
#'
#' \dontrun{
#' A <- add1CA(6,2504,2)
#' dim(A)
#' eCAN(6,2504,2) ## much better
#' Ns(4,2504,2) ## its so much worse because this ingredient is poor
#' 3047 + (1394 - 272)*2
#' eCAN(4,2504,2) ## Colbourn and Zhou power construction not implemented
#' }
#'

#' @export
add1CA <- function(t, k, v, ...){
   N <- N_add1CA(t, k, v)
   if (!is.na(N))
   add1(bestCA(t,k-1,v),t,v)
}

#' @export
add1 <- function(D, t, v, ...){
   Call <- sys.call()
   stopifnot(is.matrix(D))
   levs <- levels.no(D)
   stopifnot(all(levs==v))
   stopifnot(is.numeric(t))
   stopifnot(t%%1==0)
   stopifnot(t>=2)
   k <- ncol(D)
   aus <- cbind(D, D[,k])
   if (t==2)
     aus <- rbind(aus,
         cbind(matrix(NA,nrow=v*(v-1),ncol=k-1),
         OD2(v)[,1:2]))
   if (t==3){
     aus <- rbind(aus,
         cbind(kronecker(matrix(1,nrow=v*(v-1),ncol=1),matrix(0:(v-1),v,k-1)),
         kronecker(OD2(v)[,1:2],matrix(1,nrow=v,ncol=1))))
   }
   if (t>3){
     aus <- rbind(aus,
         cbind(kronecker(matrix(1,nrow=v*(v-1),ncol=1),bestCA(t-2,k-1,v)),
         kronecker(OD2(v)[,1:2],matrix(1,bestN(t-2,k-1,v),ncol=1))))
   }
   attrs <- attributes(D)
   return(aus)
}

#' @export
N_add1CA <- function(t, k, v){
   hilf <- bestN(t, k-1, v)
   #if (bestN(t,k,v, exclude="add1CA") <= hilf) return(NA)
   if (t==2) return(hilf + v*(v-1))
   if (t==3) return(hilf + v^2*(v-1))
   ## now t>3
   hilf + unname(bestN(t-2, k-1, v))*v*(v-1)
}
#'

## unfinished and very experimental
add2 <- function(D, t, v, ...){
  Call <- sys.call()
  stopifnot(is.matrix(D))
  levs <- levels.no(D)
  stopifnot(all(levs==v))
  stopifnot(is.numeric(t))
  stopifnot(t%%1==0)
  stopifnot(t>=2)
  k <- ncol(D)
  aus <- cbind(D, D[,k], D[,k])
  print("hier0")
  if (t==2){
    hilf <- maxconstant(bestCA(2,3,v), remove=TRUE)
    aus <- rbind(aus,
                 cbind(matrix(NA,nrow=nrow(hilf),ncol=k-1),
                       hilf))
  }
  if (t==3){
    hilf <- maxconstant(bestCA(3,3,v), remove=TRUE)
    if (ncol(OD2(v))<3) stop("not applicable")
    aus <- rbind(aus,
                 cbind(kronecker(matrix(1,nrow=v*(v-1),ncol=1),matrix(0:(v-1),v,k-1)),
                       kronecker(OD2(v)[,1:3],matrix(1,nrow=v,ncol=1))),
                 cbind(matrix(NA,nrow=nrow(hilf),ncol=k-1),
                       hilf))
  }
  if (t>3){
    if (ncol(OD2(v))<3) stop("not applicable")
    ## t-tuples with triples of the last three columns
    hilf2 <- maxconstant(bestCA(3,3,v), remove=TRUE)
    hilf1 <- bestCA(t-3,k-1,v)
print("hier")
    ## t-tuples with pairs of the last three columns
    hilf22 <- maxconstant(bestCA(2,3,v), remove=TRUE)
    hilf21 <- bestCA(t-2,k-1,v)
print("hier2")
    aus <- rbind(aus,
                 cbind(kronecker(matrix(1,nrow=nrow(hilf22),ncol=1),hilf21),
                       kronecker(hilf22, matrix(1,nrow=nrow(hilf21),ncol=1))),
                 cbind(kronecker(matrix(1,nrow=nrow(hilf2),ncol=1),hilf1),
                       kronecker(hilf2, matrix(1,nrow=nrow(hilf1),ncol=1)))
                 )
  }
  attrs <- attributes(D)
  attrs$dim[2] <- attrs$dim[2]+2
  return(aus)
}

