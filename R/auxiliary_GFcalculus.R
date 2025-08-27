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

mygf <- function(q, primitive){
  if (!q %in% 8:9) stop("mygf is only for 8 and 9")
  ## the value of primitive matters as follows,
  ##   q=8: 13 is the default, anything else is the alternative
  ##   q=9: 14 is the default, 10 and 17 are two different alternative,
  ##      of which 17 is needed; i.e., for the moment,
  ##        14 is the default, anything else yields the alternative 17

  gf <- lhs::create_galois_field(q)
  defaultprimitive <- q+sum(gf_minus(rep(0,gf$n),gf$xton, gf)*gf$p^(0:(gf$n-1)))
  ## 13 for 8, 14 for 9
  if (primitive==defaultprimitive){
    mygf <- gf
  }else{
    ## the addition and multiplication tables have been typed in from the
    ## Mathematica demonstrator
    mygf <- gf
    if (q==8){
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
        mygf$times <- rbind(rep(0,9),
                        0:8,
                        c(0,2,1,6,8,7,3,5,4),
                        c(0,3,6,4,7,1,8,2,5),
                        c(0,4,8,7,2,3,5,6,1),
                        c(0,5,7,1,3,8,2,4,6),
                        c(0,6,3,8,5,2,4,1,7),
                        c(0,7,5,2,6,4,1,8,3),
                        c(0,8,4,5,1,6,7,3,2))

      mygf$neg <- sapply(0:8, function(obj) which(mygf$plus[obj+1,]==0) - 1)
      mygf$inv <- sapply(1:8, function(obj) which(mygf$times[obj+1,]==1) - 1)
      mygf$xton <- NULL
      mygf$root <- NULL
    }
  }
  mygf
}
