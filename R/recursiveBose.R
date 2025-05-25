#' function to create a strength 2 CA based on recursive multiplication of Bose CAs
#'
#' recBoseCA applies productCA or productPCA recursively to Bose CAs. The number of levels must be a prime or prime power.
#'
#' @rdname recursiveBose
#'
#' @aliases recBoseCA
#' @aliases recursiveBose
#' @aliases N_k_recursiveBose
#' @aliases N_d_recursiveBose
#' @aliases Dd
#' @aliases kd
#'
#' @usage recBoseCA(t=2, k=NULL, v=4, type="PCA", ...)
#' @usage recursiveBose(q, d=NULL, k=NULL, type="PCA", ...)
#' @usage N_k_recursiveBose(q, d=2, type="PCA", ...)
#' @usage N_d_recursiveBose(q, k, type="PCA", ...)
#' @usage Dd(d, q)
#' @usage kd(d, q)
#'
#' @param t the requested strength; must be 2
#' @param k \code{NULL}, or the number of columns;\cr
#'          if \code{NULL}, the maximum possible number for \code{d} is returned
#'          (with \code{d} defaulting to 2, if not specified).
#' @param v the same as \code{q}
#' @param q a prime or prime power number of levels
#' @param type \code{"CA"} or \code{"PCA"},
#'          indicating whether \code{\link{productCA}} or
#'          \code{\link{productPCA}} should be used for the product
#' @param d \code{NULL}, or a small positive integer;\cr
#'          \code{d=1} returns \code{SCA_Bose(q)} (q^2 x (q+1)),
#'          \code{d} > 1 the product (\code{\link{productCA}} or \code{\link{productPCA}}, depending on \code{type})
#'          of \code{d} \code{SCA_Bose(q)};\cr
#'          For \code{d=NULL}, the minimum \code{d} for accommodating \code{k} is used.\cr
#'          If \code{k} is also \code{NULL}, the default is \code{d=2}.
#' @param ... currently not used
#'
#' @section Details:
#' The function \code{recBoseCA} uses function \code{recursiveBose}, which in turn uses
#' \code{\link{SCA_Bose}} for creating the \code{q^2 x (q+1)} ingredient matrix,
#' and either function \code{\link{productPCA}} or \code{\link{productCA}}
#' for obtaining the recursive product (depending on the choice for argument \code{type}).
#'
#' \code{type=PCA}: For a given \code{v}, the function can construct
#' up to \code{kd(d,v)} columns in \code{d*v^2-(d-1)*v} runs,
#' where \code{d} is the number of Bose CAs that are combined,
#' and \code{k(d,v)} has been determined by experimentation; for \code{d=2,3},
#' \code{k(d,v)} coincides with the value stated in Colbourn et al. (2006) Lemma 3.5,
#' which uses a recursively defined \code{Dd(0,v)=0}, \code{Dd(1,v)=1}, \code{Dd(d+1,v)=v*(Dd(d,v) + Dd(d-1,v))}
#' to yield \code{k=(v+1)*Dd(d,v) + v*Dd(d-1,v)}. Different from the construction
#' underlying Lemma 3.5 of Colbourn et al. (2006), this function improves the
#' size of the left-hand side partition in each interim step; this is achieved by
#' applying function \code{\link{CA_to_PCA}} after moving the first block of \code{v} rows
#' Experimentation showed that this choice yielded the overall optimum result among
#' blocks of \code{v} rows in all cases that were investigated.
#'
#' For \code{type=="CA"}: For a given \code{v}, the function can construct up to \code{(v+1)^d}
#' columns in up to \code{d*v^2-(d-1)*2} runs (\code{type=="CA"}),
#' where \code{d} is the number of Bose CAs that are combined.
#'
#' @returns \code{recBoseCA} returns a uniform CA strength 2 CA with \code{k} \code{v} level columns
#' (matrix of class \code{ca}); its rowsize can be determined a-priori with function \code{\link{N_recBoseCA}}.\cr
#' \code{recursiveBose} returns a uniform strength 2 CA with \code{q} level columns.\cr
#' For \code{type="PCA"} (default), the result has \code{d*q^2 - (d-1)*q} rows and
#' \code{q^d + d*q^(d-1)} columns, levels are coded from 0 to \code{q - 1}.\cr
#' For \code{type="CA"}, the result has \code{d*q^2 - (d-1)*2} rows (upper bound,
#' may be slightly smaller depending on constant rows) and
#' \code{(q+1)^d} columns, levels are coded from 0 to \code{q - 1}.\cr
#' \code{N_k_recursiveBose} returns a named vector with the number of rows and the number of columns k
#' for given \code{q} and \code{d}.\cr
#' \code{N_d_recursiveBose} returns a named vector with the run size \code{N}, the \code{d} for its creation for the requested number of columns \code{k}, and the maximum possible number of columns \code{kmax} for that construction.\cr
#'
#' @references Colbourn et al. (2006)
#'
#' @examples
#' D <- recBoseCA(k=12)   ## v=4 is the default, and type="PCA"
#' dim(D)
#' coverage(D,2)
#' eCAN(2, 12, 4)
#' Ns(2, 12, 4) ## the optimum CA is available
#'
#' D <- recBoseCA(k=25)   ## v=4 is the default, and type="PCA"
#' dim(D)
#' coverage(D,2)
#' eCAN(2, 25, 4)
#' Ns(2, 25, 4) ## the optimum CA is available
#' ## type="CA" gets close to optimal here
#' ## because k=25 makes it still work with d=2 SCA_Bose
#' D <- recBoseCA(k=25, type="CA")
#' dim(D)
#'
#' dim(D <- recursiveBose(3))
#' eCAN(2,15,3)  ## optimal
#' dim(recursiveBose(3, d=3))
#' eCAN(2,57,3)  ## reasonably close to optimal
#' dim(recursiveBose(3, d=3, type="CA"))
#' eCAN(2,64,3)  ## reasonably close to optimal
#'
#' ## query functions
#' N_k_recursiveBose(4)      ## d=2
#' N_k_recursiveBose(4, 3)   ## d=3
#' kd(3,4)   ## work horse function for k for default type PCA
#' N_k_recursiveBose(4, 3, type="CA")
#' N_d_recursiveBose(4, 12)  ## k=12
#'
#' @export
recBoseCA <- function(t=2, k=NULL, v=4, type="PCA", ...){
  Call <- sys.call()
  stopifnot(t==2)
  if (!(v %in% primedat$q)) stop("v must be prime or prime power")
  stopifnot(type %in% c("CA", "PCA"))
  suppressMessages(
    if (is.null(k))
      k <- unname(N_k_recursiveBose(v, type=type)["k"])
    else
     stopifnot(is.numeric(k), k>0, k%%1==0)
    )
  aus <- recursiveBose(v, k=k, type=type)
  class(aus) <- c("ca", class(aus))
  attr(aus, "Call") <- Call
  attr(aus, "t") <- 2
  aus
}
#'
#' @export
recursiveBose <- function(q, d=NULL, k=NULL, type="PCA", ...){
  Call <- sys.call()
  start0 <- TRUE ## default
  if (!is.null(d)){
    stopifnot(is.numeric(d))
    stopifnot(d>0)
    stopifnot(d %% 1 == 0)
  }
  if (is.null(d) && is.null(k)) d <- 2
  if (is.null(d)){
    ## infer d from k, which is not NULL
    stopifnot(is.numeric(k))
    stopifnot(k>0)
    stopifnot(k %% 1 == 0)
    if (type=="CA") d <- ceiling(log(k, base=q+1)) else{
      d <- unname(N_d_recursiveBose(q, k, type=type)["d"])
    }
  }
    if (d==1) return(SCA_Bose(q))
    A <- SCA_Bose(q)
    for (i in 2:d){
      if (type=="CA") A <- productCA(A, SCA_Bose(q)) else{
        A <- productPCA(A, SCA_Bose(q))
        nA <- nrow(A)
        # k1s <- rep(NA, nA/q)
        ## all attempts yielded the best block as the d-th block
        # for (j in 1:(nA/q)){
        #   pick <- (j-1)*q + (1:q)
        #   rest <- setdiff(1:nA, pick)
        #   hilf <- A[c(pick,rest),]
        #   k1s[j] <- attr(CA_to_PCA(hilf), "PCAstatus")$k1
        # }
        jbest <- i # which.max(k1s)
        #print(c(jbest=jbest, k1max=k1s[jbest]))
        pickbest <- (jbest-1)*q + (1:q)
        restbest <- setdiff(1:nA, pickbest)
          A <- CA_to_PCA(A[c(pickbest, restbest),])
          # A <- swapvals(A, (k1s[jbest]+1):ncol(A), 0, A[1,k1s[jbest]+1])
          # makes it an SCA; does not change anything
          # print(attr(A, "PCAstatus")$k1)
          # print(dim(A))
      }
    }
    # k1s <- rep(NA, nA/q)
    # # all attempts yielded the best block as the d-th block
    # for (j in 1:(nA/q)){
    #   pick <- (j-1)*q + (1:q)
    #   rest <- setdiff(1:nA, pick)
    #   hilf <- A[c(pick,rest),]
    #   k1s[j] <- attr(CA_to_PCA(hilf), "PCAstatus")$k1
    # }
    # jbest <- which.max(k1s)
    # pickbest <- (jbest-1)*q + (1:q)
    # restbest <- setdiff(1:nA, pickbest)
    # A <- CA_to_PCA(A[c(pickbest, restbest),])
    # print(jbest)

  class(A) <- c("ca", class(A))
  attr(A, "Call") <- Call
  attr(A, "origin") <- c("recursiveBose", type)
  attr(A, "t") <- 2
  if (!is.null(k)) A <- A[,1:k]
  A
}

#' @export
N_k_recursiveBose <- function(q, d=2, type="PCA", ...){
  if (!q %in% primedat$q)
    stop("recursiveBose requires q to be a prime or prime power")
  stopifnot(type %in% c("PCA", "CA"))
  ## check for d
  stopifnot(is.numeric(d))
  stopifnot(d > 0)
  stopifnot(d %% 1 == 0)
  if (type=="PCA"){
    return(c(N=d*q^2-(d-1)*q, k=kd(d,q)))
  } else{
    message("N can be slightly smaller because of constant rows")
    return(c(N=d*q^2-(d-1)*2, k=(q+1)^d))
  }
}

#' @export
N_d_recursiveBose <- function(q, k, type="PCA", ...){
  if (!q %in% primedat$q)
    stop("recursiveBose requires q to be a prime or prime power")
  stopifnot(type %in% c("PCA", "CA"))
  ## check for k
  stopifnot(is.numeric(k))
  stopifnot(k>1)
  stopifnot(k %% 1 == 0)
  ## figure out the necessary d
  if (type=="CA") d <- ceiling(log(k, base=q+1)) else{
    kmaxs <- sapply(1:11, function(obj) kd(obj, q))
    d <- min(which(kmaxs >= k))
    kmax=kmaxs[d]
  }
  if (type=="PCA") return(c(N=d*q^2-(d-1)*q, d=d, kmax= kmax)) else{
    message("N can be slightly smaller because of constant rows")
    return(c(N=d*q^2-(d-1)*2, d=d, kmax=(q+1)^d))
  }
}

#' @export
Dd <- function(d, q){
  ## calculates multiplier for Lemma 3.5 of Colbourn et al. 2006
  ## d corresponds to r
  Dm2=0; Dm1=1;
  if (d==0) return(0)
  if (d==1) return(1)
  for (i in 2:d){
    cur <- q*(Dm2 + Dm1);
    Dm2 <- Dm1;
    Dm1 <- cur
  }
  return(cur)
}

#' @export
kd <- function(d,q){
  ## correct split for implementation of recursiveBose
  ## with last q rows for P-block for creation
  ## and subsequent improvement at each interim step
  ## same split as for Colbourn et al. (2006) Lemma 3.5
  ## for d=2,3: (q+1)*Dd(d,q) + q*Dd(d-1,q)
  ## for d+1: k calculable from d as q*k_d + k1_d
  ##
  if (!q %in% primedat$q) stop("q must be prime or prime power")
  if (d<=3){
    k1 <- (q+1)*Dd(d,q)
    return(k1 + q*Dd(d-1,q))
  }else{
    if (d==4) return(q*((q+1)*Dd(3,q) + q*Dd(3-1,q)) + (q+1)*Dd(3,q))
    else{
      if (q > 10) return(NA)  ## huge cases (all k at least 20000 for d<=4) not implemented
      if (q==2) return(switch(as.character(d), NA,
                              '5'=164, '6'=448, '7'=1224, '8'=3344, '9'=9136,
                              '10'=24960))
      if (q==3) return(switch(as.character(d), NA,
                              '5'=828, '6'=3168, '7'=11988, '8'=45900))
      if (q==4) return(switch(as.character(d), NA,
                              '5'=2756, '6'=13280, '7'=64400))
      if (q==5) return(switch(as.character(d), NA,
                              '5'=7105, '6'=42005))
      if (q==7) return(ifelse(d==5, 31073, NA))
      if (q==8) return(ifelse(d==5, 56584, NA))
      if (q==9) return(ifelse(d==5, 96561, NA))
    }
  }
}

## test cases for d= 2, 3, 4, >4, >10


