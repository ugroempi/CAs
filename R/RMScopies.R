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
#' @aliases RMScopies_tested
#'
#' @usage RMScopiesCA(k, v, ...)
#' @usage RMScopies(q, m=0, tested_only=TRUE)
#' @usage N_RMScopiesCA(t=3, k, v, tested_only=TRUE, ...)
#' @usage RMScopies_tested
#'
#' @param k integer: the target number of columns
#' @param v integer: the number of levels (must be a prime or prime power)
#' @param q integer: a prime or prime power
#' @param m integer: the number of added columns (shifted copies),
#'          \code{0 <= m <= q^2+q+1} (\code{m=0} for \code{q=2})
#' @param t integer: the target strength (always 3)
#' @param tested_only logical: if \code{TRUE} (default), only the combinations
#'          \code{(q, m)} listed in \code{RMScopies_tested} are allowed.
#'          \code{RMScopies(q, m, tested_only=FALSE)} builds other combinations
#'          for exploration (e.g. as examples or counterexamples for the proof);
#'          such arrays should be checked with \code{caverify::ca_verify}.
#' @param ... currently not used
#'
#' @details
#' \bold{Tested parameters}: the data frame \code{RMScopies_tested} lists all
#' combinations \code{(q, m)} for which the array produced by this R code was
#' verified exhaustively (all triples of columns, all value combinations).
#' \code{RMScopiesCA(k, v)} and \code{N_RMScopiesCA} only use these combinations:
#' for k columns they take the smallest tested number of copies that is at least
#' \code{k - (v^2+v+1)}, and keep
#' the first k columns (dropping columns of a CA yields a CA). Extend the data
#' frame after verifying further cases.
#' The proof that the construction is valid for all (q, m) has not yet been
#' checked independently; until it is, only the tested combinations are used
#' by default.
#'
#' All arithmetic is in the Galois field GF(q), as provided by
#' \code{\link[lhs]{create_galois_field}} (elements coded as 0, ..., q-1;
#' for prime q, this is arithmetic modulo q).
#'
#' \bold{Base array} (\code{m=0}, Raaphorst, Moura and Stevens 2014):
#' the k0 = q^2+q+1 columns correspond to the points \eqn{p_j}{p_j} of the
#' projective plane PG(2,q), ordered as the powers of a root of a primitive
#' cubic polynomial over GF(q) (Singer cycle).
#' The array consists of a row of zeroes, a block B1 with entries
#' \eqn{x \cdot p_j}{x.p_j} for all non-zero x in GF(q)^3, and a block B2
#' constructed in the same way with the column order reversed.
#' Three columns that are collinear in B1 are not collinear in B2, so that
#' all triples of columns are covered.
#'
#' \bold{Shifted copies} (\code{m>=1}): new column \eqn{n_i}{n_i} is a copy of
#' column \eqn{c_i = i}{c_i = i} in the zero row and in B1, and equals
#' \eqn{c_i + 1}{c_i + 1} (in GF(q)) in B2. The only uncovered 3-tuples are then
#' found in column triples \eqn{\{c_i, n_i, w\}}{{c_i, n_i, w}}, for differences
#' \eqn{n_i - c_i}{n_i - c_i} different from 0 and 1 (field elements coded 2, ..., q-1).
#'
#' \bold{Layers} (\code{m>=2}): each layer consists of q^2 runs indexed by
#' y in GF(q)^2. The q+1 directions of PG(1,q) yield q+1 mutually orthogonal
#' linear functions \eqn{y \cdot P}{y.P}. In a layer, up to q copies are active:
#' active copy i receives its own direction (\eqn{c_i = y \cdot P_s}{c_i = y.P_s},
#' \eqn{n_i = y \cdot P_s + e}{n_i = y.P_s + e}),
#' all other columns receive direction \eqn{P_0}{P_0}. Each copy is active in
#' q-2 layers, with e running through the field elements coded 2, ..., q-1. A final run (0 for the original columns,
#' 1 for the new columns) covers the image of x=0 in B2.
#' The number of layers is \eqn{L = \max(q-2, \lceil m(q-2)/q \rceil)}{L = max(q-2, ceiling(m(q-2)/q))}.
#'
#' For \code{m=1}, the q(q-2) layers' worth of runs are replaced by q^2(q-2) simple runs
#' (column \eqn{c_1=a}{c_1=a}, column \eqn{n_1 = a+e}{n_1 = a+e}, all other columns equal to s),
#' and no final run is needed.
#'
#' Run sizes: \eqn{2q^3-1}{2q^3-1} for \code{m=0},
#' \eqn{2q^3-1+q^2(q-2)}{2q^3-1+q^2(q-2)} for \code{m=1},
#' \eqn{2q^3+q^2 L}{2q^3+q^2 L} for \code{m>=2}.
#'
#' Compared to the Colbourn tables (status November 2024), this improves
#' 78 entries: v=5 (k=32,...,37), v=7 (k=58,...,65), v=8 (k=74,...,81),
#' v=9 (k=92,...,101), v=11 (k=134,...,145), v=13 (k=184,...,197),
#' v=19 (k=382,...,401); for example CA(324; 3, 32, 5) instead of 349 runs,
#' CA(350; 3, 37, 5) instead of 365 runs, and CA(2025; 3, 100, 9) instead of 2169 runs.
#'
#' @returns \code{RMScopies_tested} is a data frame with columns \code{q},
#' \code{m}, \code{k} (= q^2+q+1+m), \code{N} (run size) and \code{verifier}.\cr
#' Function \code{RMScopiesCA} returns a matrix of class \code{ca}
#' with attributes, which is a strength 3 covering array.\cr
#' Function \code{RMScopies} returns a matrix without class and attributes,\cr
#' function \code{N_RMScopiesCA} returns the run size of the CA
#' (\code{NA} if the construction is not applicable).
#'
#' @references Raaphorst, S., Moura, L. and Stevens, B. (2014). A construction for
#' strength-3 covering arrays from linear feedback shift register sequences.
#' Designs, Codes and Cryptography 73, 949-968.
#'
#' @examples
#' D <- RMScopies(5)        ## the RMS array CA(249; 3, 31, 5)
#' dim(D)
#' E <- RMScopiesCA(32, 5)  ## CA(324; 3, 32, 5)
#' dim(E)
#' N_RMScopiesCA(3, 36, 5)  ## 325
#' N_RMScopiesCA(3, 37, 5)  ## 350
#' head(RMScopies_tested)
#' ## an untested combination is refused unless explicitly requested
#' try(RMScopies(11, 50))
#' dim(RMScopies(11, 50, tested_only=FALSE))
#' F <- RMScopiesCA(22, 4)  ## prime power: CA(159; 3, 22, 4)
#' dim(F)
#' \dontrun{
#' caverify::ca_verify(E, 3)$covered
#' }
#'
#' @export
RMScopies <- function(q, m=0, tested_only=TRUE){
  stopifnot(length(q) == 1, q == round(q), q >= 2)
  if (!.RMS_isprimepower(q)) stop("q must be a prime or prime power")
  k0 <- q^2 + q + 1
  stopifnot(m == round(m), m >= 0, m <= k0)
  if (q == 2 && m > 0) stop("for q=2, only m=0 is supported")
  if (tested_only && !.RMS_tested(q, m))
    stop("(q, m) = (", q, ", ", m, ") is not in RMScopies_tested; ",
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
  if (!.RMS_isprimepower(v)) stop("v must be a prime or prime power")
  k0 <- v^2 + v + 1
  m <- max(0, k - k0)
  if (m > k0) stop("k must not exceed 2*(v^2+v+1)")
  m <- .RMS_tested_m(v, m)   ## smallest tested m' >= m (columns are then dropped)
  if (is.na(m))
    stop("no tested RMScopies array for k = ", k, ", v = ", v,
         " (see RMScopies_tested)")
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
  if (!.RMS_isprimepower(v)) return(NA)
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

## internal helpers
.RMS_N <- function(q, m){
  ## run size of RMScopies(q, m)
  if (m == 0) return(2*q^3 - 1)
  if (m == 1) return(2*q^3 - 1 + q^2*(q - 2))
  2*q^3 + q^2*max(q - 2, ceiling(m*(q - 2)/q))
}

.RMS_tested <- function(q, m)
  any(RMScopies_tested$q == q & RMScopies_tested$m == m)

.RMS_tested_m <- function(q, m){
  ## smallest tested number of copies m' >= m for this q (NA if none)
  cand <- RMScopies_tested$m[RMScopies_tested$q == q & RMScopies_tested$m >= m]
  if (length(cand) == 0) NA else min(cand)
}

.RMS_isprimepower <- function(n){
  if (length(n) != 1 || !is.numeric(n) || n != round(n) || n < 2) return(FALSE)
  p <- 2
  while (n %% p != 0) p <- p + 1     ## smallest prime factor
  while (n %% p == 0) n <- n / p
  n == 1
}

## element-wise sum / product in GF(q), keeping the shape of x
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

#' @rdname RMScopies
#' @export
RMScopies_tested <- local({
  ## Combinations (q, m) whose arrays, as produced by RMScopies(q, m), were
  ## verified exhaustively (October 2026):
  ## - all except q = 19: the R output was checked with ca_verify.c
  ##   (C kernel of caverify, from github.com/chrispbsmolen/ca-tools);
  ## - q = 19: the R output coincides (same multiset of rows) with the arrays
  ##   of the Python implementation, which passed an independent brute-force
  ##   checker (verify.c in github.com/MichelBacon/covering-arrays-rms-extension).
  aus <- rbind(
    data.frame(q = 3,  m = 0:13),
    data.frame(q = 4,  m = 0:21),
    data.frame(q = 5,  m = 0:31),
    data.frame(q = 7,  m = 0:57),
    data.frame(q = 8,  m = 0:73),
    data.frame(q = 9,  m = 0:91),
    data.frame(q = 11, m = c(0:13, 22)),
    data.frame(q = 13, m = 0:14),
    data.frame(q = 16, m = c(0, 1, 16)),
    data.frame(q = 19, m = c(1, 19, 20))
  )
  aus$k <- aus$q^2 + aus$q + 1 + aus$m
  aus$N <- mapply(.RMS_N, aus$q, aus$m)
  aus$verifier <- ifelse(aus$q == 19, "identical to verified Python array",
                         "ca_verify (ca-tools)")
  aus
})
