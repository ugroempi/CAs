#' Functions to construct a strength 4 or strength 5 CA from several input CAs
#'
#' based on Martirosyan and Trung (2004). The functions obtain a strength 4
#' or strength 5 CA (strength 5 not yet working)
#' in 2k columns at v levels each from several smaller arrays.
#'
#' @rdname Martirosyan
#'
#' @aliases Martirosyan4
#' @aliases Martirosyan5
#'
#' @usage Martirosyan4(A, B, C, D, start0=TRUE, verify=FALSE, ...)
#' @usage Martirosyan5(A, B, C, D, start0=TRUE, verify=FALSE, ...)
#'
#' @param A a strength 4 or 5 CA in k columns at v levels (0, ..., v-1 or 1, ..., v);
#' for k<=3, A can be smaller (see Details section).
#'
#' @param B a strength 3 or 4 CA in k columns at v levels (0, ..., v-1 or 1, ..., v)
#'
#' @param C a strength 2 or 3 CA in k columns at v levels (0, ..., v-1 or 1, ..., v);
#' constant rows can be omitted from C, and at least one rows can always be made
#' constant by permuting levels of some columns.
#'
#' @param D a strength 2 CA in v columns at v levels (0, ..., v-1 or 1, ..., v).
#' For v=2, D can be chosen smaller (see Details section).
#'
#' @param start0 logical: do the values start at 0? This must hold for all input arrays.
#'
#' @param verify logical: should the properties of input matrices be verified?
#' (this may take a long time for larger arrays)
#'
#' @param ... further arguments to function \code{\link{coverage}}, in case
#' \code{verify} is \code{TRUE} (esp., \code{parallel argument})
#'
#' @returns a strength 4 CA with 2k columns at v levels each (same coding as ingoing CAs)
#'        in \code{N_A + (v-1)*N_B + 2*N_D*N_C} rows, where N_array denotes the run
#'        sizes of the arrays
#'
#' @details
#' For \code{Martirosyan4}: The matrices need the lower of the two strengths, i.e., 4,3,2,2 for A,B,C,D.
#' C and D can in some cases be chosen smaller than CAs:
#' One can omit constant rows from \code{C} (and one can recode it to make at least one row constant),
#' and for v=2, one can choose matrix
#' \code{D} as 2x2 with constant rows (one row for each level). For k<=3, \code{A} can be chosen identical to \code{B}.
#'
#' For \code{Martirosyan5}: The matrices need the higher of the two strengths, i.e., 5,4,3,2 for A,B,C,D.
#'
#' Unless \code{verify} is TRUE, the user is responsible for providing suitable matrices. If the matrices are not adequate, the resulting array may be worse than expected.
#'
#' @section Warning: At present, \code{Martirosyan5} does not yield strength 5, but only almost strength 5,
#' the root cause has not yet been found. Martirosyan and Trung (2004) have a typo in function h, which should
#' have an obvious correction (replace first or with and). Further mistakes in the paper or the code
#' have so far not been found.
#'
#'
#' @examples
#'   A <- lhs::createBusht(5,6,4, bRandom=FALSE)
#'   B <- lhs::createBusht(5,6,3, bRandom=FALSE)
#'   C <- lhs::createBusht(5,6,2, bRandom=FALSE)
#'   D <- C[,-6]  ## need only v=5 columns
#'   C <- C[-1,]  ## constant first row of C omitted
#'   E <- Martirosyan4(A, B, C, D)
#'   coverage(E, 4)
#'   dim(E)
#'   eCAN(4, 12, 5)  ## E is far from optimal
#'
#'   A <- lhs::createBusht(4,5,5, bRandom=FALSE)
#'   B <- lhs::createBusht(4,5,4, bRandom=FALSE)
#'   C <- lhs::createBusht(4,5,3, bRandom=FALSE)
#'   D <- lhs::createBusht(4,4,2, bRandom=FALSE)
#'   E <- Martirosyan5(A, B, C, D)
#'   coverage(E, 5)  ## something wrong
#'                   ## mistake in paper?
#'                   ## bug in programming?
#'   dim(E)
#'   eCAN(4, 12, 5)  ## run size of E is far from optimal,
#'                   ## even though not strength 5
#'

#' @export
Martirosyan4 <- function(A, B, C, D, start0=TRUE, verify=FALSE, ...){
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

#' @export
Martirosyan5 <- function(A, B, C, D, start0=TRUE, verify=FALSE, ...){
  message("Martirosyan5 does not yield strength 5, only almost strength 5. Root cause to be found.")
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



Martirosyant <- function(..., start0=TRUE, verify=FALSE){
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

