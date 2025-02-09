#' Paley constructions
#'
#' Function to construct CAs in 2-level columns with
#'     number of runs a multiple of 4,
#'     number of columns one less than runs,
#'     and coverage strength 3. The arrays are also OAs of strength 2.
#'     Also function for creating a conference matrix,
#'     and overview function for
#'
#' @aliases paley
#' @aliases show_paley
#' @aliases conference_paley
#'
#' @usage paley(q)
#' @usage show_paley(minN=NULL, maxN=NULL)
#' @usage conference_paley(q)
#'
#' @param q prime or prime power
#' @param minN,maxN bounds for the run size N. Without specified bounds,
#'        the returned data frame has more than 1000 rows.
#'
#' @section Details:
#' Based on finite fields and conference matrices,
#' the Paley construction yields arrays in numbers
#' of levels a multiple of 4. Constructions depend on
#' the value of \code{q %% 4}, which can either be 1 or 3:\cr
#' For \code{q %% 4 = 3}, a (q+1) x q CA is obtained,\cr
#' \code{q %% 4 = 1} yields a (2q+2) x (2q+1) CA.\cr
#' It may be of interest that the arrays are also orthogonal
#' arrays of strength 2.
#'
#' \code{show_paley} shows the run sizes that can be obtained
#' with the constructions, and the corresponding values of k and q, as well as
#' the smallest number for a strength 3 array in the same number of
#' 2-level columns according to the Colbourn tables
#' (currently-known upper bound for CAN) with the corresponding source entry,
#' and the run size of the NIST catalogued array.
#'
#' @returns Function \code{conference_paley} returns a conference matrix,
#' function \code{paley} a strength 3 CA (matrix) in \code{q} or \code{2q+1} 2-level columns
#' (see Section Details). Function \code{show_paley} returns a data frame
#' with columns N, k, q, constr, CAN, Source (source entry in the Colbourn tables
#' for the design with fewest possible runs), and N_NISTcat for q up to 10009.
#' Most of the arrays listed can be created;
#' exceptions are q that are powers of primes larger than 50.
#'
#' @examples
#' # paley construction I for N=12
#' D12I <- paley(11)
#' # paley construction II for N=12
#' D12II <- paley(5)
#' # they are not the same
#' table(D12I - D12II)
#' # but they have the some coverage properties
#' coverage(D12I, 3)
#' coverage(D12I, 4)
#' coverage(D12II, 3)
#' coverage(D12II, 4)
#' # paley(13) and paley(27) have different coverages
#' # for t=4
#'
#' show_paley(maxN=36)
#'   # there is no paley construction for N=16
#'   # but there is a CA(3,14,2) in 16 runs based on a
#'   # Hadamard matrix in 16 runs, which can be obtained
#'   # from package FrF2 as pb(16, nfactors=14)
#' show_paley(minN=100, maxN=150)
#'

#' @export
#' @rdname paley
conference_paley <- function(q){
  # Ball, W. W. R. and Coxeter, H. S. M.
  # Mathematical Recreations and Essays,
  # 13th ed. New York: Dover, pp. 308-309, 1987.
  gf <- lhs::create_galois_field(q)
  squares <- sort(unique(diag(gf$times)))[-1]
  nonsquares <- setdiff(1:(q-1),squares)
  x <- rep(0,q)
  x[squares+1] <- 1
  x[nonsquares+1] <- -1

  C <- matrix(0,q,q)  ## jakobsthal
  for (a in 0:(q-1))
    for (b in 0:(q-1))
      C[a+1,b+1] <- x[gf_minus(a,b,gf)+1]
  C
}

#' @export
paley <- function(q){
  stopifnot(q>=5)
  Q <- conference_paley(q)
  ## Paley construction I (q%%4==3)
  ## with 0 and 1 (q+1)x(q+1)
  ## 7, 11, 19, 23, 27, 31, 43, 47, 59, 67, 71, 79, 83, 103, 107, 127, 131, 139, 151, 163, 167, 179, 191, 199, 211, 223, 227, 239, 251, 263, 271, 283, 307, 311, 331, 347, 359, 367, 379, 383
  ## 8, 12, 20, 24, 28, 32, 44, 48, 60, 68, 72, 80, 84, 104, 108, 128, 132, 140, 152, 164, 168, 180, 192, 200, 212, 224, 228, 240, 252, 264, 272, 284, 308, 312, 332, 348, 360, 368, 380, 384
  if (q %% 4 ==3){
    hilf <- rbind(c(0,rep(1,q)),cbind(-1,Q))
    H <- t(diag(q+1)+hilf)[,-1]
    H[H==-1]<-0
  }

  ## Paley construction II (q%%4==1)
  ## with 0 and 1  (2q+2)x(2q+2)
  ##  5, 9,  13, 17, 25, 29, 37,  41, 49,   53, 61,   73, 81, 89, 97, 101, 109, 113, 121, 137, 149, 157, 169, 173, 181, 193, 197
  ## 12, 20, 28, 36, 52, 60, 76, 84, 100, 108, 124, 148, 164, 180, 196, 204, 220, 228, 244, 276, 300, 316, 340, 348, 364, 388, 396
  if (q %% 4 ==1){
    hilf <- rbind(c(0,rep(1,q)),cbind(1,Q))
    H <- matrix(0,2*(q+1), 2*(q+1))
    for (i in 1:(q+1))
      for (j in 1:(q+1)){
        if (hilf[i,j]==0) H[(2*i-1):(2*i),(2*j-1):(2*j)] <- matrix(c(1,0,0,0),2,2) else{
          if (hilf[i,j]==1) H[(2*i-1):(2*i),(2*j-1):(2*j)] <- matrix(c(1,1,1,0),2,2) else
            H[(2*i-1):(2*i),(2*j-1):(2*j)] <- matrix(c(0,0,0,1),2,2) }
      }
    H[2,] <- 1-H[2,]
    H <- H[,-1]
  }
  H
}

#' @export
show_paley <- function(minN=NULL, maxN=NULL){
  qs_paleyI <- setdiff(primedat$q[primedat$q %% 4 == 3],c(3,7))
  ## 3 and 7 do not yield coverage strength 3, do all others? tried 31, 127, both OK
  qs_paleyII <- primedat$q[primedat$q %% 4 == 1]
  N_paleyI <- qs_paleyI+1
  N_paleyII <- 2*qs_paleyII+2
  N <- c(N_paleyI, N_paleyII)
  both <- data.frame(N=N,
                     k=N-1,
                     q=c(qs_paleyI, qs_paleyII),
                     constr=rep(c("I","II"),
                                times=c(length(qs_paleyI), length(qs_paleyII))))
  both <- both[ord(both),]
  rownames(both) <- 1:nrow(both)
  if (is.null(minN) && is.null(maxN))
    return(cbind(both, do.call("rbind", lapply(both$k, function(obj) eCAN(3, obj, 2))),
           N_NISTcat=sapply(both$k, function(obj) N_NISTcat(3, obj, 2))))
  if (!is.null(minN))
    both <- both[both$N>=minN,]
  if (!is.null(maxN))
    both <- both[both$N<=maxN,]
  cbind(both, do.call("rbind", lapply(both$k, function(obj) eCAN(3, obj, 2))),
        N_NISTcat=sapply(both$k, function(obj) N_NISTcat(3, obj, 2)))
}
