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
#' DCA: DCA according to Shokri and Moura (2025)
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
    hilf <- ffCA(3, q+1)
    if (q>1) hilf <- hilf[which(lengths(apply(hilf,1,unique))==3),] else
      hilf <- hilf[which(lengths(apply(hilf,1,unique))==2),]
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

circularDCA <- function(vec){
  ## matrix of cyclic permutations
  ## used in DCA
  k <- length(vec)
  aus <- matrix(NA,k,k)
  hilf <- vec
  aus[1,] <- hilf
  if (k>1)
    for (i in 2:k){
      hilf <- c(hilf[2:k], hilf[1])
      aus[i,] <- hilf
    }
  aus
}


### the ingredients below are for
### CCL which has not yet been included
### and may be not as good as thought

N_DCA <- function(q, k){
  if (q==6 & k<=4) return(7)
  if (q==10 & k<=5) return(11)
  if (q==18 & k<=4) return(19)
  if (!q %in% primedat$q) return(NA)
  ## prime powers
  if (k<=q) return(q)
  ## odd prime powers and k=q+1
  if (q%%2==1 && k==q+1) return(q + (q-1)/2)
  ## all other cases
  return(ceiling(log(k,q))*(q-1)+1)
}

DCA <- function(q,k){
  ## q a prime or prime power
  ## k a number of columns
  ## if k<=q use Table 3
  ## if q odd and k=q+1, use construction of ???
  ## else use direct product construction
  stopifnot(q %in% primedat$q || (q==6 && k<=4) || (q==10 & k<=5))
  if (q==6){
    aus <- cbind(
      c(0, 0, 0, 0, 0, 0, 0),
      c(0, 0, 1, 2, 3, 4, 5),
      c(0, 1, 5, 4, 3, 3, 2),
      c(0, 4, 3, 5, 2, 5, 1))
    return(aus[,1:k])
  }
  if (q==10){
    aus <- cbind(
      c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(0, 0, 9, 7, 5, 6, 8, 4, 3, 1, 2),
      c(0, 1, 9, 2, 3, 5, 4, 7, 5, 8, 6),
      c(0, 2, 4, 0, 9, 5, 6, 1, 8, 7, 3),
      c(0, 6, 4, 9, 5, 3, 7, 8, 1, 2, 5)
    )
    return(aus[,1:k])
  }
  if (q==18){
    aus <- cbind(
      c(0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,  0, 0, 0, 0,  0, 0, 0),
      c(0, 0, 1, 2, 3, 4, 5, 6, 7,  8, 9,10, 11,12,13,14, 15, 16, 17),
      c(0, 9, 2, 1, 5, 7, 3, 10,14,16, 0,15, 17, 4, 6, 8, 11, 13, 12),
      c(0, 9,12,14, 6,11, 1, 7, 13, 5, 9, 2, 10,16, 8, 4, 17, 3, 15)
    )
    return(aus[,1:k])
  }
  e <- primedat$primitive[which(primedat$q==q)]
  gf <- lhs::create_galois_field(q)
  # aus <- rbind(0,cbind(0,
  #               matrix(gf_pow(rep(e,(q-1)^2),
  #                   circularDCA(0:(q-2)), gf),
  #                   nrow=q-1))
  #              )
  aus <- gf$times
  if (k<=q) return(aus[,1:k])
  if (k==q+1 && !q%%2==0){
    Dstar <- aus[-1,-1]
    D1 <- Dstar[1:((q-1)/2),,drop=FALSE]
    D2 <- Dstar[((q-1)/2+1):(q-1),,drop=FALSE]
    # X1 <- gf_pow(rep(e, (q-1)/2), 1:((q-1)/2), gf)
    X1 <- Dstar[1:((q-1)/2),1]
    aus <- rbind(
      0,
      cbind(D1,X1,0),
      cbind(D1,0,X1),
      cbind(D2,0,0)
    )
    colnames(aus) <- NULL
    return(aus)
  }
  if ((k>q+1 || q%%2==0)){
    howmany <- ceiling(log(k,base=q))
    if (howmany==2)
      return(rbind(0,productCA_raw(aus[-1,,drop=FALSE],
                                   aus[-1,,drop=FALSE]))[,1:k])
    else hilf <- aus[-1,,drop=FALSE]
    for (j in 2:howmany){
      hilf <- productCA_raw(hilf, aus[-1,,drop=FALSE])
    }
    return(rbind(0,hilf)[,1:k])
  }
}

# sizes <- matrix(NA,8,9)
# dimnames(sizes) <- list(k=3:10, q=2:10)
# for (k in 3:10){
#   for (q in 2:10){
#     aus <- try(DCA(q,k), silent=TRUE)
#     if (!"try-error" %in% class(aus))
#       if (ncol(aus)==k) sizes[k-2,q-1] <- nrow(aus)
#   }
# }
# hilf <- sizes[,c(1:4,6:8)]
# hilf <- matrix(paste0(hilf,",",hilf-1),nrow=8)
# dimnames(hilf) <- list(k=format(3:10), q=c(2:5,7:9))
# print(hilf, quote=FALSE)
#
# DCA(4,6)
