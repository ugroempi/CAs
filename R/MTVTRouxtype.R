#' Functions to construct a strength 4 or strength 5 CA from several input CAs
#'
#' based on Section 4 of Martirosyan and van Trung (2004). The functions obtain a strength 4
#' or strength 5 CA in 2k columns at v levels each from several smaller arrays.
#'
#' @rdname MTVTRouxtype
#'
#' @aliases MTVTRouxtypeCA
#' @aliases N_upper_MTVTRouxtypeCA
#'
#' @usage MTVTRouxtypeCA(t, k, v, method=NULL, ...)
#' @usage N_upper_MTVTRouxtypeCA(t, k, v, method=NULL, theoretical=FALSE, ...)
#'
#' @param t 4 or 5: strength of the CA
#' @param k positive integer: number of columns in the result
#' @param v positive integer: number of levels
#' @param method \code{NULL}, or "large"
#' @param theoretical logical: if TRUE, N is calculated from the best known CAs, not from the best implemented ones
#' @param ... currently not used
#'
#' @returns Function \code{MTVTRouxtypeCA} returns a strength \code{t} CA with \code{k} columns at \code{v} levels each, as requested.\cr
#'        \code{N_upper_MTVTRouxtypeCA} returns the size of the constructed CA before removing duplicates.
#'
#' @section Details:
#'
#' The constructions are based on a strength \code{t} CA A, a strength {t-1} CA B, a strength \code{t-2} CA C
#' and, for strength 4 with Theorem 4.3 as well as 5, a strength 2 CA D.
#' These ingredient CAs have c columns at \code{v}, where \code{c=ceiling(k/2)}, except
#' for Theorem 4.3, where D has \code{v} columns at v levels.\cr
#'
#' The construction for strength 4 is based on Theorem 4.13 (default) or on request on Theorem 4.3 (may be better for large \code{k}).
#'
#' The construction for strength 5 is based on Theorem 4.13 (default) or on request on Theorem 4.11 (may be better for large \code{k}).
#'
#' @section Warning:
#'
#' Function \code{MTVTRouxtypeCA} implements the constructions of Section 4 from Martirosyan and van Trung (2004).
#' It uses internal workhorse functions that correspond to the different theorems of the section.\cr
#' Note that the constructions are of the so-called Roux type, i.e., recursive with several
#' CA ingredients. Such constructions are competitive for very large scenarii only, and when
#' they get good ingredients.\cr
#' WARNING: Many ingredients for large situations in this package are (at present) quite poor,\cr
#' and it is easy to bring the computer to the brink of crash or beyond by naively requesting
#' a CA of strength 5 for 2000 columns. If in doubt,
#' check the expected size with \code{N_upper_MTVTRouxtypeCA}; on the author's Windows machine,
#' the creation of a CA(172792, 5, 1000, 4) by the command \code{MTVTRouxtypeCA(5, 1000, 4)}
#' was successful after 10 minutes (about half of which were used for the creation of the
#' largest ingredient with the command \code{powerCA(5, 500, 4)}) and produced an object
#' of size 1.4GB, with large RAM usage during the slow array generation.\cr
#' Furthermore, note that the creation of many or even most large arrays requires
#' an internet connection, as arrays from the NIST (2008) library are involved as
#' ingredients (in the above case, as ingredients to the \code{powerCA} construction;
#' without an internet connection, the ingredient size would more than double).\cr
#' Large arrays are not a priority of this package, so that this situation may remain
#' unchanged for a long time.
#'
#' @references Martirosyan and van Trung (2004), NIST (2008)
#'
#' @examples
#'   dim(D <- MTVTRouxtypeCA(4, 12, 3))
#'   N_upper_MTVTRouxtypeCA(4,12,3) ## two duplicates were removed
#'   N_upper_MTVTRouxtypeCA(4,12,3, theoretical=TRUE) ## ingredients are optimal
#'   coverage(D, 4)
#'   eCAN(t=4, k=12, v=3)  ## D is far from optimal
#'
#'   ## the other method is worse for small k
#'   dim(Dworse <- MTVTRouxtypeCA(4,12,3, method="large"))
#'   N_upper_MTVTRouxtypeCA(4,12,3, method="large") ## 32 duplicates were removed
#'   coverage(Dworse, 4)
#'
#'   dim(D5 <- MTVTRouxtypeCA(5, 10, 4))
#'   N_upper_MTVTRouxtypeCA(5, 10, 4) ## 48 duplicates were removed
#'   coverage(D5, 5)
#'   eCAN(t=5, k=10, v=4)  ## run size of D5 is far from optimal
#'
#'   dim(Dlarge <- MTVTRouxtypeCA(4, 360, 3))
#'   N_upper_MTVTRouxtypeCA(4, 360, 3) ## no duplicated rows removed
#'   dim(Dlarge <- MTVTRouxtypeCA(4, 360, 3, method="large")) ## better
#'   N_upper_MTVTRouxtypeCA(4, 360, 3, method="large") ## 42 duplicated rows removed
#'   Ns(4,360,3)  ## not competitive
#'
#' \dontrun{
#' ## needs an internet connection for downloading two arrays
#' ## and also takes some time for its construction
#'   N_upper_MTVTRouxtypeCA(5, 100, 4, method="large")
#'   N_upper_MTVTRouxtypeCA(5, 100, 4) ## better
#'   N_upper_MTVTRouxtypeCA(5, 100, 4, method="large", theoretical=TRUE)
#'   N_upper_MTVTRouxtypeCA(5, 100, 4, theoretical=TRUE)
#'   Ns(5, 100, 4) ## much worse than eCAN, but better than other implemented constructions
#'   dim(DLarge <- MTVTRouxtypeCA(5, 100, 4, method="large"))
#' }
#'

#' @export
MTVTRouxtypeCA <- function(t, k, v, method=NULL, ...){
  Call <- sys.call()
  stopifnot(t %in% 4:5)
  if (!is.null(method)) stopifnot(method=="large")
  kwork <- ceiling(k/2)
  tA <- min(t, kwork)
  tB <- min(t-1, kwork)
  tC <- min(t-2, kwork)
  A <- bestCA(tA, kwork, v)
  gc()
  B <- bestCA(tB, kwork, v)
  gc()
  C <- bestCA(tC, kwork, v)
  gc()
  if (t==4){
    if (is.null(method))
      aus <- MTVTRouxtype4_Thm413(A,B,C)
    else
      aus <- MTVTRouxtype4_Thm43(A,B,C)
  }else{
    ## now t=5; only theorem 4.13 implemented
    if (is.null(method)){
      D <- bestCA(2, kwork, v)
      gc()
      aus <- MTVTRouxtype5_Thm413(A,B,C,D)
    }else{
      D <- bestCA(2, v, v)
      aus <- MTVTRouxtype5_Thm411(A,B,C,D)
    }
  }
  aus <- aus[,1:k]
  class(aus) <- c("ca", class(aus))
  attr(aus, "t") <- t

  attr(aus, "Call") <- Call
  attr(aus, "origin") <- paste0("Martirosyan and van Trung (2004) Roux type, ",
                                ifelse(is.null(method), "Thm 4.13", ifelse(t==4, "Thm 4.3", "Thm 4.11")))
  aus
}

#' @export
N_upper_MTVTRouxtypeCA <- function(t, k, v, method=NULL, theoretical=FALSE, ...){
  if (!t %in% 4:5) return(NA)
  if (!is.null(method)) if (!method=="large") stop("method can be NULL or 'large'")
  stopifnot(is.numeric(k), is.numeric(v))
  stopifnot(k%%1==0, v%%1==0)
  kwork <- ceiling(k/2)
  if (is.null(method)){
    ## Theorem 4.13
      if (t==4 && !theoretical) return(unname(bestN(4, kwork, v) +
                                            (v-1)*bestN(3, kwork, v) +
                                            bestN(2, kwork, v)^2))
      if (t==4 && theoretical) return(unname(eCAN(4, kwork, v)$CAN +
                                    (v-1)*eCAN(3, kwork, v)$CAN +
                                      eCAN(2, kwork, v)$CAN^2))
    if (t==5 && !theoretical) return(unname(bestN(5, kwork, v) + (v-1)*bestN(4, kwork, v) +
                                              2*bestN(3, kwork, v)*bestN(2, kwork, v)))
    if (t==5 && theoretical) return(unname(eCAN(5, kwork, v)$CAN +
                                             (v-1)*eCAN(4, kwork, v)$CAN +
                                             2*eCAN(3, kwork, v)$CAN*eCAN(2, kwork, v)$CAN))
  }else{
    ## Theorem 4.3
    if (t==4 && !theoretical) return(unname(bestN(4, kwork, v) +
                                              (v-1)*bestN(3, kwork, v) +
                                              2*bestN(2, v, v)*bestN(2, kwork, v)))
    if (t==4 && theoretical) return(unname(eCAN(4, kwork, v)$CAN +
                                             (v-1)*eCAN(3, kwork, v)$CAN +
                                             2*eCAN(2, v, v)$CAN*eCAN(2, kwork, v)$CAN))
    ## Theorem 4.11
    if (t==5 && !theoretical) return(unname(bestN(5, kwork, v) + (v-1)*bestN(4, kwork, v) +
                                              (6*v*(v-1) + 2*bestN(2, v, v))*bestN(3, kwork, v)))
    if (t==5 && theoretical) return(unname(eCAN(5, kwork, v)$CAN +
                                             (v-1)*eCAN(4, kwork, v)$CAN +
                                             (6*v*(v-1) + 2*eCAN(2, v, v)$CAN)*eCAN(3, kwork, v)$CAN))
  }
}


MTVTRouxtype5_Thm413 <- function(A, B, C, D, start0=TRUE, ...){
  ## based on Theorem 4.13
  ## based on applying Martiroysian's t variant to strength 5
  ## the specific strength 5 variant appears to be wrong (and uses more runs)
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
  ## eliminate a constant 0-row from C
  ## cannot use remove=TRUE, because that removes all constant rows, if there are already more
  C <- maxconstant(C, one_is_enough = TRUE)
  C <- C[-1,]

  ## FD is the family of functions f1 ... fi ... fND that picks the
  ##     ith element from the column specified by the function's argument
  E1 <- cbind(A,A)
  ## B with its cyclic permutations
  E2 <- cbind(do.call(rbind, lapply(1:(v-1), function(obj) B)),
              do.call(rbind, lapply(1:(v-1), function(obj) (B-1+obj)%%v + 1)))
  ## combine C and D as indicated in Martirosyan
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
  aus <- aus[!duplicated(aus),]
  aus
}

MTVTRouxtype4_Thm413 <- function(A, B, C, start0=TRUE, ...){
  ## based on Theorem 4.13
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
  ## combine C with itself (i = t-i for t=4 and i=2) as indicated in Martirosyan
  ## i=2: slowchanging C at the top left, slowchanging C at the bottom left
  ##      fastchanging C at the top right, fastchanging C at the bottom right
  E3 <- cbind(C[rep(1:nrow(C),each=nrow(C)),],
                    C[rep(1:nrow(C),nrow(C)),])
  aus <- rbind(E1,E2,E3,cbind(E3[,(k+1):(2*k)], E3[,1:k]))
  if (start0) aus <- aus-1
  aus <- aus[!duplicated(aus),]
  aus
}

MTVTRouxtypet <- function(..., start0=TRUE){
  ### one needs a flexible list of matrices as arguments,
  ### At,...,A2
  ### no high priority - do different things first
  ## not (yet?) implemented

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
  ## check which As
  E1 <- cbind(A5,A5)
  E2 <- cbind(do.call(rbind, lapply(1:(v-1), function(obj) A4)),
              do.call(rbind, lapply(1:(v-1), function(obj) (A4-1+obj)%%v + 1)))
  E3 <- cbind(do.call(rbind, lapply(1:nrow(A4), function(obj) A3)),
              do.call(rbind, lapply(1:nrow(A4), function(obj)
                matrix(A4[obj,][A3], nrow=nrow(A3)))))
  E4 <- cbind(E3[,(k+1):(2*k)],E3[,1:k])
  aus <- rbind(E1,E2,E3,E4)
  if (start0) aus <- aus-1
  aus
}

MTVTRouxtype4_Thm43 <- function(A, B, C, start0=TRUE, ...){
  ## theorem 4.3
  ## including proposition 4.5 (which makes some matrices smaller
  ##                            for some cases)
  ## run size better for k > 2*(v+1) (v prime power)
  if (start0){
    A <- A+1
    B <- B+1
    C <- C+1
  }
  v <- max(A)
  k <- ncol(A)
  stopifnot(max(B)==v, max(C)==v)
  stopifnot(ncol(B)==k, ncol(C)==k)

  D <- bestCA(2, v, v) + 1   ## strength 2 CA
  ## proposition 4.5
    ## proposition 4.5 claims that C does not need its constant rows
    ## however, coverage is destroyed if they are removed
    #C <- maxconstant(C, t=2, remove=TRUE)
    ## for k < 4, A=B is sufficient
    if (k< 4) A <- B
    ## for 2 levels, D can be 2x2
    if (v==2) D <- matrix(1:2,2,2)

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
  unique(aus)  ## removes duplicated rows
}

MTVTRouxtype5_Thm411 <- function(A, B, C, D, start0=TRUE, ...){
  if (start0){
    A <- A+1
    B <- B+1
    C <- C+1
    D <- D+1
  }
  v <- max(A)
  if (v<3) stop("This construction requires v >= 3.")
  k <- ncol(A)
  stopifnot(max(B)==v, max(C)==v, max(D)==v)
  stopifnot(ncol(B)==k, ncol(C)==k, ncol(D)==v)
  ## FD is the family of functions f1 ... fi ... fND that picks the
  ##     ith element from the column specified by the function's argument
  E1 <- cbind(A,A) ## bestN(5, k, v)
  ## B with its cyclic permutations
  E2 <- cbind(do.call(rbind, lapply(1:(v-1), function(obj) B)),
              do.call(rbind, lapply(1:(v-1), function(obj) (B-1+obj)%%v + 1)))
        ## (v-1)*bestN(4, k, v)
  ## functions f applied, like in strength 4 construction
  E3 <- cbind(do.call(rbind, lapply(1:nrow(D), function(obj) C)),
              do.call(rbind, lapply(1:nrow(D), function(obj)
                matrix(D[obj,][C], nrow=nrow(C)))))
  E3 <- rbind(E3,
              cbind(E3[,(k+1):(2*k)],E3[,1:k]))
        ## 2*bestN(2,v,v)*bestN(3,k,v)
  ## prepare for functions g, gbar and h
  hilf <- cbind(nchoosek(v,2), nchoosek(v,2)[2:1,])
  npairs <- 2*choose(v,2)
  g <- function(x,a,b) ifelse(x==a, a, b)
  E4 <- cbind(do.call(rbind, lapply(1:npairs, function(obj) C)),
              do.call(rbind, lapply(1:npairs, function(obj){
                a <- hilf[1,obj]; b <- hilf[2,obj]
                g(C,a,b)
              })))

  E4 <- rbind(E4,
              cbind(E4[,(k+1):(2*k)],E4[,1:k]))
      ## 2*v*(v-1)*bestN(3, k, v)
  gquer <- function(x,a,b) ifelse(x==b, a, b)
  E5 <- cbind(do.call(rbind, lapply(1:npairs, function(obj) C)),
              do.call(rbind, lapply(1:npairs, function(obj){
                a <- hilf[1,obj]; b <- hilf[2,obj]
                gquer(C,a,b)
              })))
  E5 <- rbind(E5,
              cbind(E5[,(k+1):(2*k)],E5[,1:k]))
      ## 2*v*(v-1)*bestN(3, k, v)
  h <- function(x,a,b) ifelse(x==a | x==b, b, a)
    ## h replaces a with b and everything else except b with a
  E6 <- cbind(do.call(rbind, lapply(1:npairs, function(obj) C)),
              do.call(rbind, lapply(1:npairs, function(obj){
                a <- hilf[1,obj]; b <- hilf[2,obj]
                h(C,a,b)
              })))
  E6 <- rbind(E6,
              cbind(E6[,(k+1):(2*k)],E6[,1:k]))
      ## 2*v*(v-1)*bestN(3, k, v)
  aus <- rbind(E1,E2,E3,E4,E5,E6)

  if (start0) aus <- aus-1
  aus <- aus[!duplicated(aus),]
  aus
}
