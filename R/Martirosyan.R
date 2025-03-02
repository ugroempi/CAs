#' Functions to construct a strength 4 or strength 5 CA from several input CAs
#'
#' based on Martirosyan and Trung (2004). The functions obtain a strength 4
#' or strength 5 CA in 2k columns at v levels each from several smaller arrays.
#'
#' @rdname Martirosyan
#'
#' @aliases Martirosyan4
#' @aliases Martirosyan5
#'
#' @usage Martirosyan4(A, B, C, start0=TRUE, ...)
#' @usage Martirosyan5(A, B, C, D, start0=TRUE, ...)
#'
#' @param A a strength 4 or 5 CA in k columns at v levels (0, ..., v-1 or 1, ..., v);
#' for k<=3, A can be smaller (see Details section).
#'
#' @param B a strength 3 or 4 CA in k columns at v levels (0, ..., v-1 or 1, ..., v)
#'
#' @param C a strength 2 or 3 CA in k columns at v levels (0, ..., v-1 or 1, ..., v);
#' for strength 5, it appears that constant rows can be omitted from C, and at least one rows can always be made
#' constant by permuting levels of some columns.
#'
#' @param D a strength 2 CA in k columns at v levels (0, ..., v-1 or 1, ..., v).
#'
#' @param start0 logical: do the values start at 0? This must hold for all input arrays.
#'
#' @param ... currently not used
#'
#' @returns a strength 4 CA with 2k columns at v levels each (same coding as ingoing CAs)
#'        in \code{N_A + (v-1)*N_B + 2*N_C^2} rows,\cr
#'        or a strength 5 CA with 2k columns at v levels each in
#'        \code{N_A + (v-1)*N_B + 2*N_C*N_D} rows,\cr
#'        where N_array denotes the run
#'        sizes of the arrays
#'
#' @details
#' The functions implement the construction listed for general strength t CAs
#' in Martirosyan and Trung (2004).
#' The specifically-stated construction for t=4 also worked and yielded the same run size
#' for the examples considered.
#' The specific construction for t=5 yielded a larger run size and nevertheless did
#' not yield 5-coverage, i.e., there seems to be something wrong with it.
#'
#' For \code{Martirosyan4}: The matrices need the lower of the two strengths, i.e., 4,3 and 2 for A,B,C.
#'
#' For \code{Martirosyan5}: The matrices need the higher of the two strengths, i.e., 5,4,3,2 for A,B,C,D.
#'
#' The user is responsible for providing suitable matrices. If the matrices are not adequate, the resulting array may be worse than expected.
#'
#' @examples
#'   A <- lhs::createBusht(5,6,4, bRandom=FALSE)
#'   B <- lhs::createBusht(5,6,3, bRandom=FALSE)
#'   C <- lhs::createBusht(5,6,2, bRandom=FALSE)
#'   E <- Martirosyan4(A, B, C)
#'   coverage(E, 4)
#'   dim(E)
#'   eCAN(4, 12, 5)  ## E is far from optimal
#'
#'   A <- lhs::createBusht(4,5,5, bRandom=FALSE)
#'   B <- lhs::createBusht(4,5,4, bRandom=FALSE)
#'   C <- lhs::createBusht(4,5,3, bRandom=FALSE)
#'   D <- lhs::createBusht(4,5,2, bRandom=FALSE)
#'   E <- Martirosyan5(A, B, C[-1,], D)  ## omit constant row
#'   coverage(E, 5)
#'   dim(E)
#'   eCAN(4, 12, 5)  ## run size of E is far from optimal,
#'


#' @export
Martirosyan5 <- function(A, B, C, D, start0=TRUE, ...){
  ## fixed by applying the t variant to strength 5
  ## the other variant appears to be wrong (and uses more runs)
  if (start0){
    A <- A+1
    B <- B+1
    C <- C+1
    D <- D+1
  }
  v <- max(A)
  k <- ncol(A)
  stopifnot(max(B)==v, max(C)==v, max(D)==v)
  stopifnot(ncol(B)==k, ncol(C)==k, ncol(D)==k)  ## here: k columns needed
  ## FD is the family of functions f1 ... fi ... fND that picks the
  ##     ith element from the column specified by the function's argument
  E1 <- cbind(A,A)
  ## B with its cyclic permutations
  E2 <- cbind(do.call(rbind, lapply(1:(v-1), function(obj) B)),
              do.call(rbind, lapply(1:(v-1), function(obj) (B-1+obj)%%v + 1)))
  ## combine C and D as indicated in Martisoryan
  ## i=3: slowchanging C on the left, fastchanging D on the right
  E3 <- cbind(C[rep(1:nrow(C),each=nrow(D)),],
                    D[rep(1:nrow(D),nrow(C)),])
  ## i=2: slowchanging D on the left, fastchanging C on the right
  ## basically only swapping columns, speed of change placing
  ##       should be irrelevant
  E4 <- cbind(D[rep(1:nrow(D),each=nrow(C)),],
                    C[rep(1:nrow(C),nrow(D)),])
  aus <- rbind(E1,E2,E3,E4)
  if (start0) aus <- aus-1
  aus
}

#' @export
Martirosyan4 <- function(A, B, C, start0=TRUE, ...){
  if (start0){
    A <- A+1
    B <- B+1
    C <- C+1
  }
  v <- max(A)
  k <- ncol(A)
  stopifnot(max(B)==v, max(C)==v)
  stopifnot(ncol(B)==k, ncol(C)==k)  ## here: k columns needed
  E1 <- cbind(A,A)
  ## B with its cyclic permutations
  E2 <- cbind(do.call(rbind, lapply(1:(v-1), function(obj) B)),
              do.call(rbind, lapply(1:(v-1), function(obj) (B-1+obj)%%v + 1)))
  ## combine C with itself (i = t-i for t=4 and i=2) as indicated in Martisoryan
  ## i=2: slowchanging C at the top left, slowchanging C at the bottom left
  ##      fastchanging C at the top right, fastchanging C at the bottom right
  E3 <- cbind(C[rep(1:nrow(C),each=nrow(C)),],
                    C[rep(1:nrow(C),nrow(C)),])
  aus <- rbind(E1,E2,E3,cbind(E3[,(k+1):(2*k)], E3[,1:k]))
  if (start0) aus <- aus-1
  aus
}


Martirosyant <- function(..., start0=TRUE){
  ### one needs a flexible list of matrices as arguments,
  ### At,...,A2

  ## not yet implemented

  if (start0){
    A5 <- A5+1
    A4 <- A4+1
    A3 <- A3+1
    A2 <- A2+1
  }
  v <- max(A5)
  k <- ncol(A5)
  stopifnot(max(A4)==v, max(A3)==v, max(A2)==v)
  stopifnot(ncol(A4)==k, ncol(A3)==k, ncol(A2)==v)
  ## FD is the family of functions f1 ... fi ... fND that picks the
  ##     ith element from the column specified by the function's argument
  E1 <- cbind(A5,A5)
  E2 <- cbind(do.call(rbind, lapply(1:(v-1), function(obj) B)),
              do.call(rbind, lapply(1:(v-1), function(obj) (B-1+obj)%%v + 1)))
  E3 <- cbind(do.call(rbind, lapply(1:nrow(D), function(obj) C)),
              do.call(rbind, lapply(1:nrow(D), function(obj)
                matrix(D[obj,][C], nrow=nrow(C)))))
  E4 <- cbind(E3[,(k+1):(2*k)],E3[,1:k])
  aus <- rbind(E1,E2,E3,E4)
  if (start0) aus <- aus-1
  aus
}

Martirosyan4_oldDifferent <- function(A, B, C, D, start0=TRUE, ...){
  ## run size was the same as now
  if (start0){
    A <- A+1
    B <- B+1
    C <- C+1
    D <- D+1
  }
  v <- max(A)
  k <- ncol(A)
  stopifnot(max(B)==v, max(C)==v, max(D)==v)
  stopifnot(ncol(B)==k, ncol(C)==k, ncol(D)==v)
  ## FD is the family of functions f1 ... fi ... fND that picks the
  ##     ith element from the column specified by the function's argument
  E1 <- cbind(A,A)
  E2 <- cbind(do.call(rbind, lapply(1:(v-1), function(obj) B)),
              do.call(rbind, lapply(1:(v-1), function(obj) (B-1+obj)%%v + 1)))
  E3 <- cbind(do.call(rbind, lapply(1:nrow(D), function(obj) C)),
              do.call(rbind, lapply(1:nrow(D), function(obj)
                matrix(D[obj,][C], nrow=nrow(C)))))
  E4 <- cbind(E3[,(k+1):(2*k)],E3[,1:k])
  aus <- rbind(E1,E2,E3,E4)
  if (start0) aus <- aus-1
  aus
}

Martirosyan5_oldwrong <- function(A, B, C, D, start0=TRUE, ...){
  message("Martirosyan5_oldwrong does not yield strength 5, only almost strength 5. Root cause to be found.")
  if (start0){
    A <- A+1
    B <- B+1
    C <- C+1
    D <- D+1
  }
  v <- max(A)
  if (v<3) stop("This construction requires v>=3.")
  k <- ncol(A)
  stopifnot(max(B)==v, max(C)==v, max(D)==v)
  stopifnot(ncol(B)==k, ncol(C)==k, ncol(D)==v)
  ## FD is the family of functions f1 ... fi ... fND that picks the
  ##     ith element from the column specified by the function's argument
  E1 <- cbind(A,A)
  ## B with its cyclic permutations
  E2 <- cbind(do.call(rbind, lapply(1:(v-1), function(obj) B)),
              do.call(rbind, lapply(1:(v-1), function(obj) (B-1+obj)%%v + 1)))
  ## functions f applied, like in strength 4 construction
  E3 <- cbind(do.call(rbind, lapply(1:nrow(D), function(obj) C)),
              do.call(rbind, lapply(1:nrow(D), function(obj)
                matrix(D[obj,][C], nrow=nrow(C)))))
  E3 <- rbind(E3,
              cbind(E3[,(k+1):(2*k)],E3[,1:k]))
  ## prepare for functions g, gbar and h
  hilf <- nchoosek(v,2)
  npairs <- choose(v,2)
  g <- function(x,a,b) ifelse(x==a, a, b)
  E4 <- cbind(do.call(rbind, lapply(1:npairs, function(obj) C)),
              do.call(rbind, lapply(1:npairs, function(obj){
                a <- hilf[1,obj]; b <- hilf[2,obj]
                g(C,a,b)
              })))

  E4 <- rbind(E4,
              cbind(E4[,(k+1):(2*k)],E4[,1:k]))
  gquer <- function(x,a,b) ifelse(x==b, a, b)
  E5 <- cbind(do.call(rbind, lapply(1:npairs, function(obj) C)),
              do.call(rbind, lapply(1:npairs, function(obj){
                a <- hilf[1,obj]; b <- hilf[2,obj]
                gquer(C,a,b)
              })))
  E5 <- rbind(E5,
              cbind(E5[,(k+1):(2*k)],E5[,1:k]))
  h <- function(x,a,b) ifelse(x==a | x==b, b, a)
  E6 <- cbind(do.call(rbind, lapply(1:npairs, function(obj) C)),
              do.call(rbind, lapply(1:npairs, function(obj){
                a <- hilf[1,obj]; b <- hilf[2,obj]
                h(C,a,b)
              })))
  E6 <- rbind(E6,
              cbind(E6[,(k+1):(2*k)],E6[,1:k]))
  aus <- rbind(E1,E2,E3,E4,E5,E6)
  if (start0) aus <- aus-1
  aus
}
