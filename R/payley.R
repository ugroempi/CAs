#' Paley constructions
#'
#' Functions to construct Paley-based CAs as well as
#' to construct further objects related to a Paley construction
#'
#' @rdname paley
#'
#' @aliases paleyCA
#' @aliases paleyHad
#' @aliases show_paleyHad
#' @aliases paley_mat
#' @aliases conference_paley
#'
#' @usage paleyCA(t, k, ...)
#' @usage paleyHad(q)
#' @usage show_paleyHad(minN=NULL, maxN=NULL)
#' @usage paley_mat(q)
#' @usage conference_paley(q)
#'
#' @param k integer number: number of columns of the CA
#' @param t integer number: strength of the CA, 3 to 6
#' @param q prime or odd prime power
#' @param minN,maxN bounds for the run size N. Without specified bounds,
#'        the returned data frame has more than 1000 rows.
#' @param ... currently not used
#'
#' @section Background on Paley-based covering arrays:
#' Function \code{paleyCA} uses numerical experiments of
#' Colbourn (2015) to provide the smallest possible strength \code{t} CA
#' for 2-level columns that can be obtained using a Paley construction.
#' Results are based on numeric searches of Paley constructions regarding
#' coverage properties of Paley matrices with or without adding a row of ones,
#' as well as doubled such matrices, where doubling a matrix M means
#' vertically stacking M and 1-M, and adding an indicator column for the
#' two halves.\cr
#' Colbourn (2015) worked out
#' the minimum numbers of replicates for each t-tuple, and where these exceed 1,
#' proposed to remove rows without deteriorating coverage. This package contains
#' a data.frame \code{PALEYcat}, which contains an executable code string
#' from processing the results of Colbourn (2015). For strength 3, only
#' constructions without doubling are included, as CAs from doubling are
#' very uncompetitive (small corrections to Colbourn's Table 1 were needed
#' for primes 7, 13 (strength 3 not possible) and 17 (no row can be removed)).
#'
#' @section Details:
#' Function \code{paleyCA} provides a covering array of strength \code{t} for 2-level columns.\cr
#' Function \code{\link{N_PALEYcat}} shows the CA run sizes that can be obtained
#' with the paley construction using function \code{paleyCA},\cr
#' function \code{\link{k_PALEYcat}} shows the number of columns that can be accommodated
#' in \code{N} runs using function \code{paleyCA}.
#'
#' Function \code{conference_paley} yields a q x q conference matrix,
#' while function \code{paley} yields a q x q Paley matrix that
#' is derived from the conference matrix by replacing the "-1" values with zeroes.
#'
#' Based on finite fields and conference matrices,
#' function \code{paleyHad} yields OAs in numbers
#' of levels a multiple of 4. Hadamard constructions depend on
#' the value of \code{q %% 4}, which can either be 1 or 3:\cr
#' For \code{q %% 4 = 3}, a (q+1) x q OA is obtained,\cr
#' \code{q %% 4 = 1} yields a (2q+2) x (2q+1) OA.
#'
#' \code{show_paleyHad} shows the Hadamard matrix run sizes that can be obtained
#' with the constructions, and the corresponding values of k and q, as well as
#' the smallest number for a strength 3 array in the same number of
#' 2-level columns according to the Colbourn tables
#' (currently-known upper bound for CAN) with the corresponding source entry.
#'
#' @returns
#' Function \code{paleyCA} returns a matrix of class \code{ca} with levels 0 and 1 in
#' \code{k} columns and \code{N_PAYLEYcat(t, k, 2)} rows.\cr
#' Function \code{conference_paley} returns a conference matrix,
#' function \code{paleyHad} a Hadamard matrix of the Paley construction (matrix)
#' in \code{q} or \code{2q+1} 2-level columns
#' (see Section Details). Function \code{show_paleyHad} returns a data frame
#' with columns N, k, q, constr, CAN, Source (source entry in the Colbourn tables
#' for the design with fewest possible runs), and N_NISTcat for q up to 10009.
#' Most of the arrays listed can be created;
#' exceptions are q that are powers of primes larger than 50.
#'
#' @references Colbourn (2015), Colbourn and Keri (2009)
#'
#' @examples
#' ## examples for paleyCA
#'
#' ## Paley-based CA that is optimal
#' D <- paleyCA(t=4, k=67)
#' dim(D)
#' # coverage(D, 4) ## commented out, as it runs too long
#' eCAN(4, 67, 2)
#' # as good as the optimum cyclotomy construction
#' # (and actually the same)
#'
#' D <- paleyCA(t=3, k=10)
#' dim(D)
#' coverage(D, 3)
#' eCAN(3, 10, 2)
#' eCAK(3, 12, 2)
#' k_PALEYcat(3, 12, 2)
#' # as good as optimum (first columns of Hadamard matrix)
#'
#' D <- paleyCA(t=5, k=300)
#' dim(D)
#' eCAN(5, 300, 2)
#' eCAK(5, 359, 2)
#' k_PALEYcat(5, 359, 2)
#' ## as good as the derived design from the Colbourn / Keri construction
#'
#' D <- paleyCA(t=5, k=400)
#' dim(D)
#' eCAN(5, 400, 2) ## current best-known
#'
#' D <- paleyCA(t=5, k=80)
#' dim(D)
#' eCAN(5, 80, 2) ## current best-known slightly better
#'
#' D <- paleyCA(t=5, k=192)
#' dim(D)
#' eCAN(5, 192, 2) ## current best-known 13 fewer runs
#'
#' D <- paleyCA(t=6, k=29)
#' dim(D)
#' N_PALEYcat(t=6, k=29)
#' eCAN(6, 29, 2) ## current best-known much better
#' k_PALEYcat(t=6, N=720) ## good for larger number of factors
#' eCAN(6, 360, 2) ## as good as the current-best
#' eCAK(6, 720, 2) ## design from the Colbourn / Keri construction
#'
#' ## other Paley-related functions
#' # conference matrix (square q x q matrix with diagonal 0)
#' conference_paley(5)
#' # paley matrix (square q x q matrix with diagonal 0)
#' paley_mat(5)
#'
#' ## Hadamard matrix based OAs
#' # paley construction I for N=12
#' D12I <- paleyHad(11)
#' # paley construction II for N=12
#' D12II <- paleyHad(5)
#' # they are not the same
#' table(D12I - D12II)
#' # but they have the same coverage properties
#' coverage(D12I, 3)
#' coverage(D12I, 4)
#' coverage(D12II, 3)
#' coverage(D12II, 4)
#' # paleyHad(13) and paleyHad(27) have different coverages
#' # for t=4
#'
#' show_paleyHad(maxN=36)
#'   # there is no paley construction for N=16
#'   # but there is a CA(3,14,2) in 16 runs based on a
#'   # Hadamard matrix in 16 runs, which can be obtained
#'   # from package FrF2 as pb(16, nfactors=14)
#' show_paleyHad(minN=100, maxN=150)
#'

#' @export
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
paleyHad <- function(q){
  stopifnot(q>=3)  ## changed 5 to 3 May 16 2025
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
show_paleyHad <- function(minN=NULL, maxN=NULL){
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
    return(cbind(both, do.call("rbind", lapply(both$k, function(obj) eCAN(3, obj, 2)))))
  if (!is.null(minN))
    both <- both[both$N>=minN,]
  if (!is.null(maxN))
    both <- both[both$N<=maxN,]
  cbind(both, do.call("rbind", lapply(both$k, function(obj) eCAN(3, obj, 2))))
}

#' @export
paley_mat <- function(q){
  ## similar to conference matrix:
  ## the -1 entries become 1, the +1 and 0 entries become 0;
  ## equal to the matrix A from the cyclotomy construction
  stopifnot(q>=3)  ## changed from 5 to 3 on May 16 2025
  gf <- lhs::create_galois_field(q)
  ## cycvec checks v and q and retrieves the primitive
  xstart <- cycvec(2,q)
  A <- matrix(NA,q,q)
  for (i in 1:q)
    for (j in 1:q){
      A[i,j] <- xstart[gf_minus(j,i, gf)+1]
    }
  A
}

#' @export
paleyCA <- function(t, k, ...){
  Call <- sys.call()
  stopifnot(is.numeric(t))
  stopifnot(is.numeric(k))
  hilf <- PALEYcat[PALEYcat[,"t"]==t & PALEYcat[,"k"]>=k,,drop=FALSE]
  if (nrow(hilf)==0) stop("impossible combination of t and k for function paleyCA")
  zeile <- which.min(hilf[,"N"])
  aus <- eval(parse(text=hilf$code[zeile]))[,1:k]
  class(aus) <- c("ca", class(aus))
  attr(aus, "Call") <- Call
  attr(aus, "origin") <- paste0("Paley construction, q=", hilf$q[zeile])
  attr(aus, "t") <- t
  aus
}
