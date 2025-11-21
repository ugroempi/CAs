#' auxiliary functions for Galois field arithmetics
#'
#' non-visible functions to be used in other functions,
#' including the construction for OA(q^3, 3, q+2, q) for q=2^s
#'
#' @rdname auxiliary_GFcalculus
#'
#' @usage NULL
#'
#' @details Most functions are taken from SOAs. Function createBush3PowerOf2
#' implements the construction via permutation vectors as stated in Sherwood (2006)
#'
#' @importFrom lhs create_galois_field

int2poly <- function(x, gf){
  gf$poly[x+1, , drop=FALSE]
}

gf_prod <- function(x, y, gf){
  #stopifnot("GaloisField" %in% class(gf))
  #q <- gf$q
  #stopifnot(all(c(x,y) %in% 0:(q-1)))
  #stopifnot(length(x)==length(y))
  sapply(seq_along(x), function(obj)
    gf$times[x[obj]+1, y[obj]+1])
}

gf_matmult <- function (M1, M2, gf, checks = TRUE)
{
  if (checks) {
    stopifnot("GaloisField" %in% class(gf))
    q <- gf$q
    stopifnot(all(c(M1, M2) %in% 0:(q - 1)))
    stopifnot(is.matrix(M1))
    stopifnot(is.matrix(M2))
    stopifnot(ncol(M1) == nrow(M2))
  }
  nc1 <- ncol(M1)
  nr1 <- nrow(M1)
  summanden <- vector(mode="list")
  for (i in 1:nc1)
    summanden[[i]] <- outer(M1[,i], M2[i,], gf_prod, gf)
  aus <- gf_sum_list(summanden, gf, checks=FALSE)
  matrix(aus, nrow=nr1)
}

gf_sum <- function(x, y, gf){
  #stopifnot("GaloisField" %in% class(gf))
  #q <- gf$q
  #stopifnot(all(c(x,y) %in% 0:(q-1)))
  #stopifnot(length(x)==length(y))
  sapply(seq_along(x), function(obj)
    gf$plus[x[obj]+1, y[obj]+1])
}

gf_minus <- function(x,y, gf){
  ## calculates x - y, after reducing both by a mod operation to
  ## gf entries
  q <- gf$q; p <- gf$p
  x <- x%%q; y <- y%%q
  yn <- gf$neg[y+1]
  mapply(gf_sum, x, yn, MoreArgs= list(gf=gf))
}

gf_sum_list <- function (ll, gf, checks = TRUE)
{
  ## ll is a list of integer vectors to be summed over gf
  if (checks) {
    stopifnot("GaloisField" %in% class(gf))
    if (!all(c(unlist(ll)) %in% 0:(gf$q - 1)))
      stop("invalid numbers occur in ll")
    if (!length(unique(lengths(ll))) == 1)
      stop("all elements of ll must have the same length")

  }
  hilf <- lapply(ll, int2poly, gf=gf)
  hilf <- base::Reduce("+", hilf)%%gf$p
  apply(hilf, 1, function(obj) lhs::poly2int(gf$p, gf$n, obj))
}

gf_prod_list <- function (ll, gf, checks = TRUE)
{
  ## ll is a list of integer vectors to be summed over gf
  if (checks) {
    stopifnot("GaloisField" %in% class(gf))
    if (!all(c(unlist(ll)) %in% 0:(gf$q - 1)))
      stop("invalid numbers occur in ll")
    if (!length(unique(lengths(ll))) == 1)
      stop("all elements of ll must have the same length")

  }
  gfprodnow <- function(obj1, obj2) CAs:::gf_prod(obj1, obj2, gf)
  base::Reduce(gfprodnow, ll)
}

gf_pow <- function(x, pow, gf, checks=TRUE){
  ## x, pow scalar
  if (checks) {
    stopifnot("GaloisField" %in% class(gf))
    if (!all(x %in% 0:(gf$q - 1)))
      stop("invalid numbers occur in x")
    if (!length(x) == length(pow))
      stop("x and pow must have the same length")
    if (any(pow<0)) stop("pow must have non-negative entries")
    if (any(!pow%%1==0)) stop("pow must have integer-valued entries")
  }
  sapply(seq_along(x), function(obj)
    ifelse(pow[obj]==0, 1, ifelse(pow[obj]==1, x[obj],
    gf_prod_list(rep(x[obj], pow[obj]), gf))))
}

gf_det <- function(mat, gf){
  d <- nrow(mat)
  stopifnot(ncol(mat)==d)
  stopifnot("GaloisField" %in% class(gf))
  stopifnot(all(mat %in% 0:(gf$q-1)))
  mat <- cbind(mat,mat[,1:(d-1)])
  ## plus-products minus minus-products
  plusprod <- gf_sum_list(lapply(0:(d-1), function(obj){
    ## obj is the column plus over the row
    gf_prod_list(sapply(1:d, function(obj2){
      ## obj2 is the row
      mat[obj2, obj+obj2]}), gf)}), gf)
  minusprod <- gf_sum_list(lapply(0:(d-1), function(obj){
    ## obj is the column plus over d minus the row + 1
    gf_prod_list(sapply(d:1, function(obj2){
      ## obj2 is the row
      mat[obj2, d-obj2+1+obj]}), gf)}), gf)
  gf_minus(plusprod, minusprod, gf)
}

createBush3PowerOf2 <- function(q){
  stopifnot(((q/2)%%1)==0)
  gf <- lhs::create_galois_field(q)
  ## N <- q^3; k <- q+2

  ## low powers at the top
  ## the h tuples for non-reduced permutation vectors
  ## q = 2^ell
  lspace <- 0:(q-1)
  lspace <- rbind(1,lspace,
                  gf_prod(lspace, lspace, gf))
  betaspace <- sfsmisc::digitsBase(0:(q^3-1), q, ndigits=3)[3:1,]
  ## the h tuples for reduced permutation vectors w.r.t. t=4
  ## (i.e., using only the rows for powers 0 to 2 from t=4 digits
  ##        or all rows from t=3 digits)
  lspaceReduced <- cbind(c(0,1,0), c(0,0,1))
  cbind(gf_matmult(t(betaspace), lspace, gf),
        gf_matmult(t(betaspace), lspaceReduced, gf))
}

mygf <- function(q, type="2018"){
  ## with version 0.22, changed type="default" to type="2018" (CL_SCPHFs)
  ##     with alternatives "2006" for SMC_SCPHFs
  ##     and               "2009" for WC_SCPHFs
  if (!q %in% c(8:9,16,25)) stop("mygf is only for 8, 9, 16, 25")
  if (type=="2009" && !q %in% 8:9) stop("type 2009 is for q=8 or 9 only")
  ## since version 0.22, there are three versions
  ## 2018: lhs::create_galois_field (as of Art Owen)
  ## 2006: SageMath / Python function numbthy of
  ##    https://github.com/Robert-Campbell-256/Number-Theory-Python
  ## 2009: https://github.com/robbywalker/ca-phf-research/blob/master/tabu/math/galois_field.cpp
  ##    this one directly gives multiplication tables (not polynomials)
  ## 2009 is for 8 and 9 only,
  ##      for 8 it is the same as 2006,
  ##      for 9 it has a different multiplication table

  ## the remark below is kept for now, but is LIKELY IRRELEVANT
  ## if more versions show up, this can be handled via a primitive
  ## the value of primitive matters as follows,
  ##   q=8: 13 is the default, anything else is the alternative
  ##   q=9: 14 is the default, 10 and 17 are two different alternatives,
  ##      of which 17 is needed; i.e., for the moment,
  ##        14 is the default, anything else yields the alternative 17
  ##   q=16 and q=25: not yet worked out, default must be worked out from xton component of
  ##                  lhs::create_galois_field outcome which gives coeffs of the rhs of x^n = a_n-1*x^(n-1) + ... + a_1*x + a_0
  ##   according to Michael Wagner, the Python version uses Conway polynomials

  gf <- lhs::create_galois_field(q)
#  defaultprimitive <- q+sum(gf_minus(rep(0,gf$n),gf$xton, gf)*gf$p^(0:(gf$n-1)))
  ## 13 for 8, 14 for 9
  if (type=="2018" || !q %in% c(8:9, 16, 25)){
    mygf <- gf
  }else{
    ## the addition and multiplication tables have been typed in from the
    ## Mathematica demonstrator or taken from robbywalker github
    mygf <- gf
    if (q==8){
      ## 2006 and 2009 use the same gf
      ## mygf$plus is unchanged vs the lhs version
      mygf$times <- rbind(
        rep(0,8),
        0:7,
        c(0, 2,	4,	6,	3,	1,	7,	5),
        c(0, 3,	6,	5,	7,	4,	1,	2),
        c(0, 4,	3,	7,	6,	2,	5,	1),
        c(0, 5,	1,	4,	2,	7,	3,	6),
        c(0, 6,	7,	1,	5,	3,	2,	4),
        c(0, 7,	5,	2,	1,	6,	4,	3)
      )
      mygf$neg <- sapply(0:7, function(obj) which(mygf$plus[obj+1,]==0) - 1)
      mygf$inv <- sapply(1:7, function(obj) which(mygf$times[obj+1,]==1) - 1)
      mygf$xton <- NULL
      mygf$root <- NULL
    }
    if (q==9){
      # default polynomial: x^2 + x + 2
      ## try x^2+1 (implies xton 2 0, i.e. primitive 10)
      ## and x^2 + 2x + 2 (implies xton 1 1, i.e., primitive 17)
      ## it seems that the addition table is the same for all polynomials
      hilf <- rbind(0:(q-1),
                    c(1,2,0,4,5,3,7,8,6),
                    c(2,0,1,5,3,4,8,6,7))
      mygf$plus <- rbind(hilf,
                    cbind(hilf[,4:9],hilf[,1:3]),
                    cbind(hilf[,7:9],hilf[,1:6]))
      ## multiplication table for polynomial x^2 + 1
      # if (primitive==10){
        # mygf$times <- rbind(rep(0,9),
        #                 0:8,
        #                 c(0,2,1,6,8,7,3,5,4),
        #                 c(0,3,6,2,5,8,1,4,7),
        #                 c(0,4,8,5,6,1,7,2,3),
        #                 c(0,5,7,8,1,3,4,6,2),
        #                 c(0,6,3,1,7,5,2,8,5),
        #                 c(0,7,5,4,2,6,8,3,1),
        #                 c(0,8,4,7,3,2,5,1,6))
      #}

      ## multiplication table for polynomial x^2 + 2x + 2
      # if (primitive==17)
      if (type=="2006")   ## for SMC
        mygf$times <- rbind(rep(0,9),
                        0:8,
                        c(0,2,1,6,8,7,3,5,4),
                        c(0,3,6,4,7,1,8,2,5),
                        c(0,4,8,7,2,3,5,6,1),
                        c(0,5,7,1,3,8,2,4,6),
                        c(0,6,3,8,5,2,4,1,7),
                        c(0,7,5,2,6,4,1,8,3),
                        c(0,8,4,5,1,6,7,3,2))
        ## for Walker and Colbourn cphfs
        ## from https://github.com/robbywalker/ca-phf-research/blob/master/tabu/math/galois_field.cpp

        if (type=="2009")
         mygf$times <- rbind(
          rep(0,9),
          0:8,
          c(0, 2, 1, 6, 8, 7, 3, 5, 4),
          c(0, 3, 6, 2, 5, 8, 1, 4, 7),
          c(0, 4, 8, 5, 6, 1, 7, 2, 3),
          c(0, 5, 7, 8, 1, 3, 4, 6, 2),
          c(0, 6, 3, 1, 7, 4, 2, 8, 5),
          c(0, 7, 5, 4, 2, 6, 8, 3, 1),
          c(0, 8, 4, 7, 3, 2, 5, 1, 6)
        )

      mygf$neg <- sapply(0:8, function(obj) which(mygf$plus[obj+1,]==0) - 1)
      mygf$inv <- sapply(1:8, function(obj) which(mygf$times[obj+1,]==1) - 1)
      mygf$xton <- NULL
      mygf$root <- NULL
    }
    if (q==16){
      ## for SMC
      mygf$times <- rbind(
      c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
      c(0, 2, 4, 6, 8, 10, 12, 14, 3, 1, 7, 5, 11, 9, 15, 13),
      c(0, 3, 6, 5, 12, 15, 10, 9, 11, 8, 13, 14, 7, 4, 1, 2),
      c(0, 4, 8, 12, 3, 7, 11, 15, 6, 2, 14, 10, 5, 1, 13, 9),
      c(0, 5, 10, 15, 7, 2, 13, 8, 14, 11, 4, 1, 9, 12, 3, 6),
      c(0, 6, 12, 10, 11, 13, 7, 1, 5, 3, 9, 15, 14, 8, 2, 4),
      c(0, 7, 14, 9, 15, 8, 1, 6, 13, 10, 3, 4, 2, 5, 12, 11),
      c(0, 8, 3, 11, 6, 14, 5, 13, 12, 4, 15, 7, 10, 2, 9, 1),
      c(0, 9, 1, 8, 2, 11, 3, 10, 4, 13, 5, 12, 6, 15, 7, 14),
      c(0, 10, 7, 13, 14, 4, 9, 3, 15, 5, 8, 2, 1, 11, 6, 12),
      c(0, 11, 5, 14, 10, 1, 15, 4, 7, 12, 2, 9, 13, 6, 8, 3),
      c(0, 12, 11, 7, 5, 9, 14, 2, 10, 6, 1, 13, 15, 3, 4, 8),
      c(0, 13, 9, 4, 1, 12, 8, 5, 2, 15, 11, 6, 3, 14, 10, 7),
      c(0, 14, 15, 1, 13, 3, 2, 12, 9, 7, 6, 8, 4, 10, 11, 5),
      c(0, 15, 13, 2, 9, 6, 4, 11, 1, 14, 12, 3, 8, 7, 5, 10)
      )
    }
    if (q==25){
      ## for SMC
      mygf$times <- rbind(
        c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24),
        c(0, 2, 4, 1, 3, 10, 12, 14, 11, 13, 20, 22, 24, 21, 23, 5, 7, 9, 6, 8, 15, 17, 19, 16, 18),
        c(0, 3, 1, 4, 2, 15, 18, 16, 19, 17, 5, 8, 6, 9, 7, 20, 23, 21, 24, 22, 10, 13, 11, 14, 12),
        c(0, 4, 3, 2, 1, 20, 24, 23, 22, 21, 15, 19, 18, 17, 16, 10, 14, 13, 12, 11, 5, 9, 8, 7, 6),
        c(0, 5, 10, 15, 20, 8, 13, 18, 23, 3, 11, 16, 21, 1, 6, 19, 24, 4, 9, 14, 22, 2, 7, 12, 17),
        c(0, 6, 12, 18, 24, 13, 19, 20, 1, 7, 21, 2, 8, 14, 15, 9, 10, 16, 22, 3, 17, 23, 4, 5, 11),
        c(0, 7, 14, 16, 23, 18, 20, 2, 9, 11, 6, 13, 15, 22, 4, 24, 1, 8, 10, 17, 12, 19, 21, 3, 5),
        c(0, 8, 11, 19, 22, 23, 1, 9, 12, 15, 16, 24, 2, 5, 13, 14, 17, 20, 3, 6, 7, 10, 18, 21, 4),
        c(0, 9, 13, 17, 21, 3, 7, 11, 15, 24, 1, 5, 14, 18, 22, 4, 8, 12, 16, 20, 2, 6, 10, 19, 23),
        c(0, 10, 20, 5, 15, 11, 21, 6, 16, 1, 22, 7, 17, 2, 12, 8, 18, 3, 13, 23, 19, 4, 14, 24, 9),
        c(0, 11, 22, 8, 19, 16, 2, 13, 24, 5, 7, 18, 4, 10, 21, 23, 9, 15, 1, 12, 14, 20, 6, 17, 3),
        c(0, 12, 24, 6, 18, 21, 8, 15, 2, 14, 17, 4, 11, 23, 5, 13, 20, 7, 19, 1, 9, 16, 3, 10, 22),
        c(0, 13, 21, 9, 17, 1, 14, 22, 5, 18, 2, 10, 23, 6, 19, 3, 11, 24, 7, 15, 4, 12, 20, 8, 16),
        c(0, 14, 23, 7, 16, 6, 15, 4, 13, 22, 12, 21, 5, 19, 3, 18, 2, 11, 20, 9, 24, 8, 17, 1, 10),
        c(0, 15, 5, 20, 10, 19, 9, 24, 14, 4, 8, 23, 13, 3, 18, 22, 12, 2, 17, 7, 11, 1, 16, 6, 21),
        c(0, 16, 7, 23, 14, 24, 10, 1, 17, 8, 18, 9, 20, 11, 2, 12, 3, 19, 5, 21, 6, 22, 13, 4, 15),
        c(0, 17, 9, 21, 13, 4, 16, 8, 20, 12, 3, 15, 7, 24, 11, 2, 19, 6, 23, 10, 1, 18, 5, 22, 14),
        c(0, 18, 6, 24, 12, 9, 22, 10, 3, 16, 13, 1, 19, 7, 20, 17, 5, 23, 11, 4, 21, 14, 2, 15, 8),
        c(0, 19, 8, 22, 11, 14, 3, 17, 6, 20, 23, 12, 1, 15, 9, 7, 21, 10, 4, 18, 16, 5, 24, 13, 2),
        c(0, 20, 15, 10, 5, 22, 17, 12, 7, 2, 19, 14, 9, 4, 24, 11, 6, 1, 21, 16, 8, 3, 23, 18, 13),
        c(0, 21, 17, 13, 9, 2, 23, 19, 10, 6, 4, 20, 16, 12, 8, 1, 22, 18, 14, 5, 3, 24, 15, 11, 7),
        c(0, 22, 19, 11, 8, 7, 4, 21, 18, 10, 14, 6, 3, 20, 17, 16, 13, 5, 2, 24, 23, 15, 12, 9, 1),
        c(0, 23, 16, 14, 7, 12, 5, 3, 21, 19, 24, 17, 10, 8, 1, 6, 4, 22, 15, 13, 18, 11, 9, 2, 20),
        c(0, 24, 18, 12, 6, 17, 11, 5, 4, 23, 9, 3, 22, 16, 10, 21, 15, 14, 8, 2, 13, 7, 1, 20, 19))
    }
  }
  mygf
}
