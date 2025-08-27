#' auxiliary ingredients
#'
#' non-visible functions for creating ingredients
#'
#' @rdname auxiliary_ingredients
#'
#' @usage NULL
#'
#' @section Details:
#'
#' OD2: ordered designs of strength 2, for the Sherwood mixed level construction
#' OD3: ordered designs of strength 3, for ordered design (CCL)
#'
## returns the bottom left part of SCA_Bose
## or a trivial OD for v not a prime power
## with exceptions for v=6 and v=10
## OD3 for strength 3 is in auxiliary_functions

OD2 <- function(v){
  ## determine the OD according to Sherwood (2008)
  ## for prime power v, it is an LOD, otherwise an OD
  ## (OD is sufficient; OD can be verified by verifying
  ## strength 2 for the OA with all-constant rows added
  ## by calculating the GWLP
  if (v %in% primedat$q) return(SCA_Bose(v)[(v+1):(v^2), -(v+1)])
  if (v==6) return(miscCA(2,3,6)[-(1:6),]) ## three columns possible
  if (v==10) return(maxconstant(miscCA(2,3,10))[-(1:10),]) ## three columns possible
  one <- rep(0:(v-1), v-1)
  cbind(one, (one + rep(1:(v-1), each=v))%%v)
}

OD3 <- function(q){
  ## CCL, with Inf replaced by q+1
  ## stopifnot(q %in% primedat$q)
  Q <- 0:(q-1)
  ## treat the case for which q is not prime or prime power
  ## (yields three columns only)
  if (!q %in% primedat$q){
    hilf <- ffCA(3,q+1)
    hilf <- hilf[which(lengths(apply(hilf,1,unique))==3),]
    return(hilf)
  }
  ## initially only prime, later use gf_mult
  gf <- lhs::create_galois_field(q)
  INV <- gf$inv[-1] # does not include position for 0
  MULT <- gf$times ## includes row and column for 0
  ADD <- gf$plus
  tuples1 <- expand.grid(1,Q,Q,Q)
  tuples1 <- tuples1[which(!tuples1[,4]==sapply(1:nrow(tuples1), function(obj) MULT[tuples1[obj,2]+1,tuples1[obj,3]+1])),]
  tuples01 <- expand.grid(0,1:(q-1),1,Q)
  tuples <- rbind(tuples1, tuples01)
  X <- matrix(NA, nrow(tuples), q+1)
  for (i in 1:nrow(X)){
    a <- tuples[i,1];
    b <- tuples[i,2];
    c <- tuples[i,3];
    d <- tuples[i,4]
    xvec <- c(Q, Inf) ## column labels
    for (j in 1:(q+1)){
      xnow <- xvec[j]
      if (xnow==Inf){
        ## j==q+1
        if (b==0) X[i,j] <- Inf else X[i,j] <- MULT[a+1,INV[b]+1]
      }else{
        if (ADD[MULT[b+1,xnow+1]+1, d+1]==0) X[i,j] <- Inf else
          X[i,j] <- MULT[ADD[MULT[a+1, xnow+1] + 1, c+1] + 1,
                         INV[ADD[MULT[b+1, xnow+1] + 1, d+1]]+1]
      }
    }
  }
  X[X==Inf] <- q
  X
}

