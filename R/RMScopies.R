#' Strength 3 CAs from the Raaphorst-Moura-Stevens construction, extended by shifted column copies
#'
#' The Raaphorst-Moura-Stevens (2014) construction yields a CA(2q^3-1; 3, q^2+q+1, q)
#' for prime powers \code{q}. Function \code{RMScopies} extends it by \code{m} additional columns,
#' each of which is a shifted copy of an existing column, plus a few layers of
#' additional runs.
#'
#' @rdname RMScopies
#'
#' @aliases RMScopiesCA
#' @aliases RMScopies
#' @aliases N_RMScopiesCA
#' @aliases k_RMScopiesCA
#'
#' @usage RMScopiesCA(k, v, ...)
#' @usage RMScopies(q, m=0, tested_only=TRUE)
#' @usage N_RMScopiesCA(t=3, k, v, tested_only=TRUE, ...)
#' @usage k_RMScopiesCA(t=3, N, v, tested_only=TRUE, ...)
#'
#' @param k integer: the target number of columns
#' @param v integer: the number of levels (must be a prime or prime power)
#' @param q integer: a prime or prime power
#' @param m integer: the number of added columns (shifted copies),
#'          \code{0 <= m <= q^2+q+1} (\code{m=0} for \code{q=2})
#' @param t integer: the target strength (always 3)
#' @param N integer: the available number of runs
#' @param tested_only logical: if \code{TRUE} (default), only the combinations
#'          \code{(q, m)} listed in \code{\link{rmsCOPIEScat}} are allowed.
#'          \code{RMScopies(q, m, tested_only=FALSE)} builds other combinations
#'          for exploration (e.g. as examples or counterexamples for the proof);
#'          such arrays should be checked with \code{caverify::ca_verify}.
#' @param ... currently not used
#'
#' @details
#' For \code{m=0}, the construction was proven by Raaphorst, Moura and Stevens (2014).
#' The extension for \code{m>0} was developed and implemented
#' in 2026 by Mickael Bergeron with Claude. As the mathematical proof
#' of validity for all pairs \code{(q, m)} with \code{m>0}
#' has not yet been checked independently by a human, the verified parameter
#' sets stored in \code{\link{rmsCOPIEScat}} are the basis of function
#' \code{RMScopiesCA}, while function \code{RMScopies} can also apply the
#' construction for non-verified pairs on request. Mathematical detail can be found
#' in the documentation for \code{\link{rmsCOPIEScat}}.
#'
#' @returns \code{N} (run size) and \code{verifier}.\cr
#' Function \code{RMScopiesCA} returns a matrix of class \code{ca}
#' with attributes, which is a strength 3 covering array.\cr
#' Function \code{RMScopies} returns a matrix without class and attributes,\cr
#' function \code{N_RMScopiesCA} returns the run size of the CA
#' (\code{NA} if the construction is not applicable),
#' function \code{k_RMScopiesCA} returns the largest number of columns
#' achievable with at most N runs (\code{NA} if none).
#'
#' @references Raaphorst, Moura and Stevens (2014)
#'
#' @examples
#' D <- RMScopies(5)        ## the RMS array CA(249; 3, 31, 5)
#' dim(D)
#' eCAN(3,31,5) ## best according to the Colbourn tables
#' E <- RMScopiesCA(32, 5)  ## CA(324; 3, 32, 5)
#' dim(E)
#' eCAN(3,32,5) ## E is 25 runs smaller
#' N_RMScopiesCA(3, 36, 5)  ## 325
#' N_RMScopiesCA(3, 37, 5)  ## 350
#' k_RMScopiesCA(3, 325, 5)  ## 36
#' head(rmsCOPIEScat)
#' ## an untested combination is refused unless explicitly requested
#' try(RMScopies(11, 50))
#' dim(RMScopies(11, 50, tested_only=FALSE))
#' F <- RMScopiesCA(22, 4)  ## prime power: CA(159; 3, 22, 4)
#' dim(F)
#' Ns(3,22,4) ## not best but among the smaller ones
#'
#' @export
RMScopies <- function(q, m=0, tested_only=TRUE){
  stopifnot(length(q) == 1, q == round(q), q >= 2)
  if (!q %in% primedat$q) stop("q must be a prime or prime power")
  k0 <- q^2 + q + 1
  stopifnot(m == round(m), m >= 0, m <= k0)
  if (q == 2 && m > 0) stop("for q=2, only m=0 is supported")
  if (tested_only && !.RMS_tested(q, m))
    stop("(q, m) = (", q, ", ", m, ") is not in rmsCOPIEScat; ",
         "use tested_only=FALSE to build it anyway, and verify the result")
  gf <- lhs::create_galois_field(q)
  ## Singer cycle: powers of a root of a primitive cubic
  P <- .RMS_singer_points(gf)[1:k0, , drop=FALSE]
  X <- as.matrix(expand.grid(0:(q-1), 0:(q-1), 0:(q-1)))[-1, , drop=FALSE]  ## non-zero vectors
  dimnames(X) <- NULL
  B1 <- .RMS_gf_mattmult(X, P, gf)
  B2 <- .RMS_gf_mattmult(X, P[k0:1, , drop=FALSE], gf)
  aus <- rbind(rep(0, k0), B1, B2)
  if (m == 0) return(aus)
  ## shifted copies of columns 1, ..., m
  inB2 <- (1 + nrow(B1) + 1):nrow(aus)
  neu <- aus[, 1:m, drop=FALSE]
  neu[inB2, ] <- .RMS_gf_plus(neu[inB2, ], 1, gf)
  aus <- cbind(aus, neu)
  k <- k0 + m
  if (m == 1){
    extra <- do.call(rbind, lapply(0:(q-1), function(a)
      do.call(rbind, lapply(2:(q-1), function(e)
        t(sapply(0:(q-1), function(s){
          r <- rep(s, k); r[1] <- a; r[k] <- .RMS_gf_plus(a, e, gf); r
        }))))))
    return(rbind(aus, extra))
  }
  ## layers
  act <- cbind(i=rep(1:m, times=q-2), e=rep(2:(q-1), each=m))
  if (m <= q) layers <- lapply(2:(q-1), function(e) cbind(i=1:m, e=e))
  else layers <- lapply(seq(1, nrow(act), by=q),
                        function(st) act[st:min(st + q - 1, nrow(act)), , drop=FALSE])
  dirs <- rbind(c(1, 0), c(0, 1), cbind(1, 1:(q-1)))   ## P_0, ..., P_q
  Y <- as.matrix(expand.grid(0:(q-1), 0:(q-1))); dimnames(Y) <- NULL
  extra <- do.call(rbind, lapply(layers, function(lay){
    stopifnot(!any(duplicated(lay[, "i"])))
    B <- matrix(.RMS_gf_mattmult(Y, dirs[1, , drop=FALSE], gf), nrow(Y), k)
    for (s in seq_len(nrow(lay))){
      f <- .RMS_gf_mattmult(Y, dirs[s + 1, , drop=FALSE], gf)[, 1]
      B[, lay[s, "i"]] <- f
      B[, k0 + lay[s, "i"]] <- .RMS_gf_plus(f, lay[s, "e"], gf)
    }
    B
  }))
  rbind(aus, extra, c(rep(0, k0), rep(1, m)))
}

#' @export
RMScopiesCA <- function(k, v, ...){
  ## t=3 always holds
  Call <- sys.call()
  if (!v %in% primedat$q) stop("v must be a prime or prime power")
  k0 <- v^2 + v + 1
  m <- max(0, k - k0)
  if (m > k0) stop("k must not exceed 2*(v^2+v+1)")
  m <- .RMS_tested_m(v, m)   ## smallest tested m' >= m (columns are then dropped)
  if (is.na(m))
    stop("no tested RMScopies array for k = ", k, ", v = ", v,
         " (see rmsCOPIEScat)")
  aus <- RMScopies(v, m)[, 1:k, drop=FALSE]
  class(aus) <- c("ca", class(aus))
  attr(aus, "origin") <- if (m == 0) "Raaphorst-Moura-Stevens 2014"
                         else "Raaphorst-Moura-Stevens 2014, with shifted column copies"
  attr(aus, "Call") <- Call
  attr(aus, "t") <- 3
  aus
}

#' @export
N_RMScopiesCA <- function(t=3, k, v, tested_only=TRUE, ...){
  if (!t == 3) return(NA)
  if (!v %in% primedat$q) return(NA)
  q <- v; k0 <- q^2 + q + 1
  m <- max(0, k - k0)
  if (m > k0) return(NA)
  if (q == 2 && m > 0) return(NA)
  if (tested_only){
    m <- .RMS_tested_m(q, m)
    if (is.na(m)) return(NA)
  }
  .RMS_N(q, m)
}

#' @export
k_RMScopiesCA <- function(t=3, N, v, tested_only=TRUE, ...){
  if (!t == 3) return(NA)
  if (!v %in% primedat$q) return(NA)
  q <- v; k0 <- q^2 + q + 1
  ms <- if (q == 2) 0 else 0:k0
  if (tested_only) ms <- intersect(ms, rmsCOPIEScat$m[rmsCOPIEScat$q == q])
  ms <- ms[sapply(ms, function(m) .RMS_N(q, m)) <= N]
  if (length(ms) == 0) return(NA)
  k0 + max(ms)
}

## internal helpers
.RMS_N <- function(q, m){
  ## run size of RMScopies(q, m)
  if (m == 0) return(2*q^3 - 1)
  if (m == 1) return(2*q^3 - 1 + q^2*(q - 2))
  2*q^3 + q^2*max(q - 2, ceiling(m*(q - 2)/q))
}

.RMS_tested <- function(q, m)
  any(rmsCOPIEScat$q == q & rmsCOPIEScat$m == m)

.RMS_tested_m <- function(q, m){
  ## smallest tested number of copies m' >= m for this q (NA if none)
  cand <- rmsCOPIEScat$m[rmsCOPIEScat$q == q & rmsCOPIEScat$m >= m]
  if (length(cand) == 0) NA else min(cand)
}

## element-wise sum / product in GF(q), keeping the shape of x;
## vectorised versions of gf_sum / gf_prod (auxiliary_GFcalculus.R),
## which loop over elements and would be slow for the large arrays (q = 19)
.RMS_gf_plus <- function(x, y, gf){
  aus <- x
  aus[] <- gf$plus[cbind(as.vector(x) + 1, as.vector(y) + 1)]
  aus
}
.RMS_gf_times <- function(x, y, gf){
  aus <- x
  aus[] <- gf$times[cbind(as.vector(x) + 1, as.vector(y) + 1)]
  aus
}

## X %*% t(P) over GF(q), for matrices X (r x d) and P (c x d)
.RMS_gf_mattmult <- function(X, P, gf){
  r <- nrow(X); c <- nrow(P)
  aus <- matrix(0L, r, c)
  for (i in seq_len(ncol(X))){
    prod <- matrix(gf$times[cbind(rep(X[, i], times=c) + 1, rep(P[, i], each=r) + 1)], r, c)
    aus <- .RMS_gf_plus(aus, prod, gf)
  }
  aus
}

.RMS_singer_points <- function(gf){
  ## returns the q^3-1 powers of alpha as rows (coefficients of 1, x, x^2),
  ## alpha a root of a primitive cubic x^3 + f2 x^2 + f1 x + f0 over GF(q)
  q <- gf$q
  mulx <- function(a, f){
    ## a: coefficients (1, x, x^2); multiply by x and reduce mod f (monic, degree 3):
    ## x^3 = -(f0 + f1 x + f2 x^2)
    .RMS_gf_plus(c(0, a[1], a[2]), gf$neg[.RMS_gf_times(rep(a[3], 3), f, gf) + 1], gf)
  }
  for (f0 in 1:(q-1)) for (f1 in 0:(q-1)) for (f2 in 0:(q-1)){
    f <- c(f0, f1, f2)
    cur <- c(1, 0, 0)
    pw <- matrix(0L, q^3 - 1, 3)
    ok <- TRUE
    for (j in 1:(q^3 - 1)){
      pw[j, ] <- cur
      cur <- mulx(cur, f)
      if (j < q^3 - 1 && all(cur == c(1, 0, 0))) { ok <- FALSE; break }
    }
    if (ok && all(cur == c(1, 0, 0)))
      return(pw)     ## alpha has order q^3-1: f is primitive
  }
  stop("no primitive polynomial found")
}
