#' Make a larger uniform CA from a DHHF and smaller CAs
#'
#' Function to construct a larger CA by replacing levels of a DHHF with columns of smaller CAs
#'
#' @rdname DHHF2CA
#'
#' @aliases DHHF2CA
#'
#' @usage DHHF2CA(P, Dlist, v = max(Dlist[[1]]) + 1, ...)
#'
#' @param P an M x k DHHF of strength t for partition size p=min(t,v) with several blocks of rows,
#'          where the ith block has w_i levels
#'          (a PHHF is the case with partition size t)
#' @param Dlist a list of N_i x w_i CAs of strength t with \code{v} levels each
#' @param v number of levels of all elements of \code{Dlist},
#'          automatically detected in case of coding as 0,...,\code{v}-1
#' @param ... further arguments to function \code{\link{maxconstant}},
#'          e.g., \code{one_is_enough=TRUE} for suppressing the search for more constant rows
#'
#' @section Details:
#' A PHHF is a perfect heterogeneous hash family, a DHHF a distributed heterogeneous hash family
#' (weaker condition); heterogeneity means that there are c blocks with different numbers of levels,
#' w_1,...,w_c, and there are u_i rows for the ith block, with u_1 +...+ u_c = M.\cr
#' Function \code{DHHF2CA} yields a CA(chi + sum_i=1^c (N_i-rho_i)*u_i, t, k, v) from
#' replacing the levels of the DHHF with the columns of the corresponding CA from
#' \code{Dlist},\cr
#' leaving out the rho_i constant rows of the ith CA and adding a few (chi) constant rows back
#' in, if needed, where chi = max(0, v − sum_i=1^c u_i*(v - rho_i)).
#' The construction is from Theorem 2.3 of Colbourn and Torres-Jimenez (2010).\cr
#' The function uses function \code{maxconstant} to maximize the constant rows of the CAs.
#'
#' There may be (likely rare) situations with many different numbers of levels in the DHHF
#' for which the function's handling of constant rows could be slightly improved;
#' for these cases, the result will not have the minimum number of rows possible from
#' the construction, but one or very few more rows.
#'
#' The function does not check for the properties of the ingoing arrays. Hence, it is the users
#' responsibility to ensure that these are as needed.
#'
#' @returns a CA(chi + sum_i=1^c (N_i-rho_i)*u_i, t, k, v)
#' (matrix, not of class \code{ca}, as there are no
#' guarantees regarding the strength, because guarding this is left to the user or to calling functions).
#'
#' @references Colbourn and Torres-Jimenez (2010), Ji and Yin (2010) for the example
#'
#' @examples
#' #########################################
#' ## DHHF2CA
#' #########################################
#'
#' ## the Ji and Yin (2010) strength 3 OA in six 12-level columns (Lemma 3.5 with Corollary 3.3)
#' ## available in this package as oa1728.12.6
#' P <- t(oa1728.12.6[,1:5])
#' P <- T(P, 4:5, 1)
#' ## two columns reduced to 11 levels each
#' dim(P)
#' eCAN(3, 12, 2) ## 15 runs, in TJ2level_CAs
#' eCAN(3, 11, 2) ## 12 runs, paleyCA also yields this
#' D <- DHHF2CA(P,
#'     list(bestCA(3,12,2), bestCA(3,11,2)))
#' dim(D)
#' eCAN(3, 1452, 2)  ## optimal, this construction
#'
#'

#' @export
DHHF2CA <- function(P, Dlist, v=max(Dlist[[1]])+1, ...){
  ## function for use by other functions
  ## skips all the checks, except for list and matrix
  stopifnot(is.list(Dlist))
  stopifnot(all(sapply(Dlist, is.matrix)))
  stopifnot(all(sapply(Dlist, is.numeric)))
  stopifnot(is.matrix(P))
  stopifnot(is.numeric(P))

  ## determine blocks of P
  mins <- apply(P,1,min)
  stopifnot(length(unique(mins))==1)
  if (mins[1]==0) P <- P+1 ## assuming no other violations
  maxs <- apply(P,1,max)
  descendingorder <- sort(maxs, decreasing = TRUE, index.return=TRUE)$ix
  P <- P[descendingorder,]
  maxs <- maxs[descendingorder]
  uis <- rev(table(maxs)) ## descending order
  wis <- as.numeric(names(uis))
  cblocks <- length(uis)

  ## check Dlist
  stopifnot(length(Dlist)==cblocks)
  wis_Dlist <- sapply(Dlist, ncol)
  descendingorder <- sort(wis_Dlist, decreasing = TRUE, index.return=TRUE)$ix
  Dlist <- Dlist[descendingorder]

  if (!(all(wis == wis_Dlist)))
    stop("The DHHF and the CAs in Dlist do not match.")

  ## which CA to pick for each row of P:
  ##

  M <- nrow(P);
  Nis <- sapply(Dlist, nrow)
  k <- ncol(P)
  levs <- sort(unique(c(Dlist[[1]])))
  if (min(levs) > 0){
    ## make levels start with 0
    Dlist <- lapply(Dlist, function(D) D - min(levs))
    levs <- levs - min(levs)
  }
  Drem <- lapply(Dlist, function(D) maxconstant(D, verbose=2, ...))
  nconst <- sapply(Drem, function(obj) length(attr(obj, "constant_rows")$row_set_list))

  ## nconst contains the rho_i of Theorem 2.3 of Colbourn and Torres-Jimenez
  ## in the order of Dlist
  chi <- max(0, v - sum(uis*(v-nconst)))
      ## if chi > 0, some constant rows have to be added back in

  const <- mapply(function(obj1, obj2) c(obj1[1:obj2,1]),
                  Drem, as.list(nconst))
       ## list of constant elements
  Drem <- mapply(function(obj1, obj2) obj1[-(1:obj2),],Drem, as.list(nconst))
       ## matrices without constant elements to use in every step

  nonconst <- lapply(const, function(obj) setdiff(levs, obj))
  aus <- matrix(NA, chi + sum(uis*(Nis-nconst)), k)
    Y <- vector(mode="list", length=M+1) ## for covered elements
    Nbisher <- 0
    Mbisher <- 0
    for (j in 1:cblocks){
      nconstnow <- nconst[j]
      constnow <- const[[j]]
      ## it might improve the situation if the matrix Drem[[j]]
      ## would be adapted such that its constant rows have as many levels
      ## as possible from coveredbisher (for j>=2)
      ## Where this is encountered, a message should complain about exceeding chi
      Ninow <- Nis[j]
      for (i in 1:uis[j]){
        aus[(Nbisher+1):(Nbisher+Ninow-nconstnow),] <-
          (Drem[[j]][,P[Mbisher+i,]] + (i-1)*nconstnow) %% v
        r <- Mbisher + i
        Y[[r]] <- setdiff(levs, (constnow + nconstnow) %% v)
        Nbisher <- Nbisher + Ninow - nconstnow
      }
      Mbisher <- Mbisher + uis[j]
      coveredbisher <- sort(unique(unlist(Y[1:Mbisher])))
  }

  ## Yr the levels that are not constant in the CA used for row r
  ## YMp1 the levels that were constant for all rows (there are such levels for chi>0)
  if (chi>0){
    YMp1 <- setdiff(levs, coveredbisher)
    if (length(YMp1) < chi) stop("wrong length for YMp1")
    if (length(YMp1) > chi){
      addrow <- length(YMp1) - chi
      warning(paste("DHHF2CA encountered a case for which the implementation is suboptimal\n",
                    addrow, " more rows than possible"))
      chi <- length(YMp1)
    aus <- rbind(aus, matrix(NA, addrow, k))
    }
    aus[(Nbisher+1):(Nbisher+chi),] <- matrix(YMp1, nrow=chi, ncol=k)
  }
  ## unique may improve the result for poorly chosen ingredients
  unique(aus)
}
