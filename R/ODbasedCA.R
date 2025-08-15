#' Ordered design based constructions by Cohen Colbourn and Ling (CCL)
#'
#' Functions to construct a strength 3 CA with k columns at v levels
#' based on the "ordered design (CCL)" construction by Cohen, Colbourn and Ling (2003)
#'
#' @aliases ODbasedCA
#' @aliases N_ODbasedCA
#' @aliases k_ODbasedCA
#'
#' @usage ODbasedCA(t, k, v, ...)
#' @usage N_ODbasedCA(t, k, v, nconstant=NULL, ...)
#' @usage k_ODbasedCA(t, N, v, ...)
#'
#' @param t integer, strength of the CA
#' @param k integer, number of columns of the CA; \code{k <= v}
#' @param v integer, number of levels of the CA; \code{v >= k}, and \code{v-1} must be a prime power
#' @param N integer, affordable run size
#' @param nconstant \code{NULL} for the default; otherwise 1 or 2, for the number of constant rows in the ingredient CA (default matches the available ingredients)
#' @param ... currently not used
#'
#' @section Details:
#' The function \code{ODbasedCA} implements constructions according to
#' Cohen, Colbourn and Ling (2003); the run size is more explicitly stated in
#' Cohen, Colbourn and Ling (2008).
#'
#' The construction works in three steps: It creates an ordered design of strength 3
#' for all triples with three distinct levels. It then uses \code{choose(v,2)} copies
#' of the best strength 3 CA in \code{k} 2-level columns (called "ingredient CA" below)
#' for covering all triples with at most two distinct levels at a time.
#'
#' The ingredient CA always has at least one constant row, and in a few cases two
#' (at present for \code{k<=8} and \code{k=12}).\cr
#' Each of its \code{choose(v,2)} copies is used for one of the \code{choose(v,2)} pairs
#' of distinct levels. If the ingredient CA has two constant rows, the run size can be optimized by
#' omitting both from all copies of the ingredient CA (i.e., by initially not covering
#' (most) triples with all elements the same) and by adding all \code{v} constant rows
#' in the end (the function actually leaves them all in and removes duplicates in the end).
#'
#' If the ingredient CA has only a single constant row, this row can be omitted entirely
#' without adding any constant rows in the end, by making sure that each of the
#' \code{v} levels occurs at least once in the role of 0 and in the role of 1 in the
#' ingredient CA (this is achieved by
#' swapping 0 and \code{v-1} in that respective pair).
#'
#' This package knows how to obtain the best available ingredient CAs (except for the missing
#' constant row for 7 and 8 2-level columns), so the run size of the output object
#' matches the Colbourn table entries (which have at least \code{k=10}). For \code{k=12},
#' the arrays are even smaller than the Colbourn table entries
#' (which do not reflect the benefit from having two constant elements in the ingredient CA).
#'
#' Function \code{N_ODbasedCA} per default calculates the run size for an ingredient
#' with two constant rows for \code{k<=6} and \code{k=12}, and with one constant row for all other cases.
#' The argument \code{nconstant} can be set to 1 or 2 in order to inspect the respective
#' non-default case.
#'
#' Function \code{k_ODbasedCA} calculates the maximum \code{k} for a given run size.
#' The upper limit for \code{k} is always \code{v}, so that this is less interesting than
#' for other constructions.
#'
#' @returns Function \code{powerCA}
#' returns a strength \code{t} CA in \code{N_powerCA(t, k, v)} runs
#' and \code{k} columns. \code{k_powerCA(t, N, v)} returns the maximum number of columns
#' that can be achieved with \code{N} runs when requesting strength \code{t} at \code{v} levels.
#'
#' @references Cohen, Colbourn and Ling (2003, 2008)
#'
#' @examples
#' D <- ODbasedCA(3, 8, 8)
#' coverage(D, 3)
#' dim(D)
#' N_ODbasedCA(3, 8, 8)
#' ## not using a 12 run design with 2 constant row
#' ## would imply 20 more rows
#' N_ODbasedCA(3, 8, 8, nconstant=1)
#' k_ODbasedCA(3, 630, 8) ## pessimistic, assumes nconstant=1
#' k_ODbasedCA(3, 650, 8) ## assumes nconstant=1
#'
#' D <- ODbasedCA(3, 12, 14)
#' coverage(D, 3)
#' dim(D)
#' N_ODbasedCA(3, 12, 14)
#' N_ODbasedCA(3, 12, 14, nconstant=1)
#' eCAN(3, 12, 14) ## based on one constant row,
#'                 ## although two are in the best-known CA(15,3,12,2)
#'

#' @export
ODbasedCA <- function(t=3, k, v, ...){
  stopifnot((v-1) %in% primedat$q)
  OD <- OD3(v-1)
  hilfN <- bestN(3, k, 2)
  ## maximum feasible k for the necessary N
  if (hilfN==12){
    hilfk <- k
    constr <- "miscCA"
  }else{
    if (names(hilfN)=="CK_doublingCA"){
      hilfk <- k
      constr <- "CK_doublingCA"
    }else{
      ## k_CK_doublingCA not in ks because of difficulties
    hilfk <- ks(3, hilfN, 2)
    hilfk <- max(hilfk[-length(hilfk)])
    constr <- names(bestN(3, hilfk, 2))
    }
  }
  ## obtain array with maximum possible number of columns
  Y <- eval(parse(text=labelToCode(constr, 3, hilfk, 2)))
  ## make at least one row constant
  Y <- maxconstant(Y)
  if (!hilfk==k){
  ## get the number of ones per row
  hilfY <- rowSums(Y)
  ## select columns with two constant rows
  ## or with as many as possible constant entries in rows
  range_hilfY <- range(hilfY)
  if (all(range_hilfY==c(0,1))) Y <- Y[,1:k] else{
    if (range_hilfY[1] > 0) Y <- 1 - Y ## now the constant row is zeroes
    ## no all-one row
    ## adjust range_hilf
    hilfY <- rowSums(Y)
    range_hilfY <- range(hilfY)
    rowpick <- which.max(hilfY)
    colpicks <- which(Y[rowpick,]==1)
    if (length(colpicks)>=k) colpicks <- colpicks[1:k] else
      colpicks <- c(colpicks, setdiff(1:hilfk, colpicks)[1:(k-length(colpicks))])
    Y <- Y[c(1,rowpick,setdiff(2:hilfN, rowpick)), colpicks]

    ## possibly achieve two constant rows
    Y <- maxconstant(Y)
  }
  }

  Y <- Y[-1,] ## remove first constant row (all zeroes)
  pairs <- nchoosek(v, 2) - 1
  ## swap, so that constant rows are all covered
  ## only violated for 0 (always first) and v-1 (always last)
  colswitch <- which(pairs[1,]==0 & pairs[2,] == v-1)
  pairs[,colswitch] <- c(v-1,0)
  ###############################################
  pairdesigns <- vector(mode="list", length=choose(v,2))
  pairdesigns[[1]] <- Y
  for (i in 2:ncol(pairs)){
    pairdesigns[[i]] <- matrix(pairs[,i][Y+1], ncol=k)
  }
  Y <- do.call(rbind, pairdesigns)
  unique(rbind(Y, OD[,1:k]))
}

#' @export
N_ODbasedCA <- function(t=3,k,v, nconstant=NULL, ...){
  stopifnot(is.numeric(t), is.numeric(k), is.numeric(v))
  stopifnot(t%%1==0, k%%1==0, v%%1==0)
  if (!is.null(nconstant)) stopifnot(nconstant %in% 1:2)
  if (!t==3) return(NA)
  if (k > v) return(NA)
  if (!(v-1) %in% primedat$q) return(NA)
  if (is.null(nconstant)) nconstant <- ifelse(k %in% c(1:8, 12), 2, 1)
  q <- v-1
  if (nconstant==2) N <- q^3-q + choose(v,2)*bestN(3, k, 2) - (q^2-1) else
    N <- q^3-q + choose(v,2)*(bestN(3, k, 2) - 1)
  unname(N)
}

#' @export
k_ODbasedCA <- function(t=3,N,v, ...){
  ## nconstant depends on k, therefore pessimistically fixed at 1
  stopifnot(is.numeric(t), is.numeric(N), is.numeric(v))
  stopifnot(t%%1==0, N%%1==0, v%%1==0)
  if (!t==3) return(NA)
  ## k<=v
  if (!(v-1) %in% primedat$q) return(NA)
  q <- v-1
  offset <- N - q^3 + q
  ## 1 constant row
    hilfN <- floor(offset/choose(v,2)+1)
    suppressMessages(hilfks <- ks(3, hilfN, 2))
    suppressWarnings(k1 <- min(v, max(hilfks[-length(hilfks)])))
  ## 2 constant rows
    hilfN <- floor((offset + q^2-1)/choose(v,2))
    suppressMessages(hilfks <- ks(3, hilfN, 2))
    suppressWarnings(k2 <- min(v, max(hilfks[-length(hilfks)])))
  if (k2 %in% c(2:6,12) && k2>k1){
    return(k2)
  }else return(ifelse(k1==-Inf, NA, k1))
}

