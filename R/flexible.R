#' Redundancy checks and profile
#'
#' Functions for manipulating flexible values: determination of flexible positions (flexpos),
#' determination of a flexible set (markflex), and calculation of the profile (flexprofile)
#'
#' @rdname flexible
#'
#' @aliases flexpos
#' @aliases markflex
#' @aliases flexprofile
#'
#' @usage flexpos(D, t, ...)
#' @usage markflex(D, t, fixrows=0, verbose=0, ...)
#' @usage flexprofile(D, ...)
#'
#' @param D an N x k CA of strength t with v levels, coded from 0 to v-1 or from 1 to v; flexible values, if any, must be denoted as \code{NA}
#' @param t integer-valued (>=2), the strength of \code{D}
#' @param fixrows integer from 0 to \code{N} for preventing the top \code{fixrows} rows from being moved
#' @param verbose integer-valued degree of verbosity; not yet implemented; intended for documenting the row permutations that were conducted
#' @param ... currently not used
#'
#' @section Details:
#' Function \code{flexpos} marks flexible \emph{positions}, i.e. positions all of whose \code{t}-way interactions are also covered
#' in at least one tuple without them. It will not work to treat all these positions as flexible \emph{at the same time},
#' because this would eliminate \emph{all} instances of the affected interactions.
#'
#' Function \code{markflex} uses the greedy algorithm of @nayeriRandomizedPostoptimizationCovering2013 for marking elements
#' as flexible. It makes use of function \code{flexpos} for swapping rows so that it becomes more likely that entire rows
#' consist of wildcards (flexible set points).
#' It works reasonably fast for strength 2 or small \code{k}, but becomes slow for larger strength and/or many columns.
#' Nayeri et al. (2013) ascertain that the determination of the maximum set of flexible points is an NP-hard problem, and
#' this algorithm is not guaranteed to find the maximum flexible set.
#' The returned flexible positions tend to be at the bottom of the array.
#'
#' Function \code{profile} obtains the profile of column-wise flexible values for \code{D}.
#'
#' @returns Function \code{flexpos} returns a logical matrix with \code{TRUE} values indicating
#' flexible positions of \code{D}. Any position that is already \code{NA} in \code{D} is set to \code{TRUE}.\cr
#' Function \code{markflex} returns the ingoing CA, with values denoted as before,
#' with rows potentially reordered (for facilitating all-flexible rows),
#' and with flexible positions changed to NA. \bold{this may change}\cr
#'
#' @references Nayeri, Colbourn and Konjevod (2013)
#'
#' @examples
#' # identify flexible values
#' A <- cyc(19,2)
#' coverage(A, 3)
#' (AwithFlex <- markflex(A, 3))
#' coverage(AwithFlex, 3)
#'
#' ## also works for mixed levels
#' (L18withFlex <- markflex(DoE.base::L18, 2))
#' flexprofile(L18withFlex)
#' ## last two rows removable
#' ## row 16 only needed for pair 23 in first two columns
#' ## change element [8,1] to 2 --> row 16 becomes removable
#' L18withFlex[8,1] <- 2
#' flexprofile(L18withFlex)  ## frequency reduced by 1 for column 1
#' (L18withFlexmod <- markflex(L18withFlex, 2))
#' ## now, three rows can be removed
#' coverage(L18withFlexmod[1:15,], 2, start0=FALSE)
#'
#' ## do the replacement on the raw design
#' L18mod <- DoE.base::L18; L18mod[8,1] <- 2
#' (L18modwithFlex <- markflex(L18mod, 2, flexible="NA"))
#'
#' ## example of Colbourn and Torres-Jimenez (2013)
#' ## Figure 2 of the paper
#' hilf <- strsplit(c("2222001", "0022222", "1121201", "0120110",
#'                         "2210202", "2102122", "1200020", "0211121",
#'                         "1012100", "2001210", "0000001", "1111012",
#'                         "2222222", "2222011"), "")
#' plan <- do.call(rbind, lapply(hilf, as.numeric))
#' dim(plan)
#' (planMod <- markflex(plan, 2))  ## last two rows can be removed
#' coverage(planMod[1:12,], 2)     ## now optimal
#' eCAN(2,7,3)
#' ## fixing the first two rows from being moved
#' ## prevents the second constant row,
#' ## i.e., fixing rows may have a price
#' (planMod2 <- markflex(plan, 2, fixrows=2))
#'
#' ## to accommodate one additional 2-level column
#'

#' @export
markflex <- function(D, t, fixrows=0, verbose=0, ...){
   hilf <- matcheck(D, uniform=FALSE, flexible=NA)
   v <- hilf$v ## vector of length k;
   k <- hilf$k; N <- hilf$N
   uniform <- FALSE
   if (length(unique(v))==1){
     uniform <- TRUE
     v <- v[1]
   }
   stopifnot(fixrows %in% 0:N)
   flexPositions <- flexpos(D, t, verbose=verbose, ...)
   ## move rows with many flexible positions to bottom
   ## goal: be able to eliminate rows
   rowOrder <- 1:N
   if (fixrows==0)
     rowOrder <- sort(rowSums(flexPositions), index.return=TRUE)$ix
   else if (fixrows < N-1)
     rowOrder <- c(1:fixrows,
                   sort(rowSums(flexPositions[(fixrows+1):N,]), index.return=TRUE)$ix + fixrows)
   D <- D[rowOrder,]

   start0 <- hilf$start0
   vstart <- as.numeric(!start0)   ## number
   vend <- v - as.numeric(start0)  ## number or vector

   flexval <- NA
   hilf <- matrix(TRUE, N, k)
   tuples <- nchoosek(k,t)
   if (uniform) {
     cands <- do.call(expand.grid, rep(list(vstart:vend), t))
     ncands <- v^t
   }
   for (i in 1:ncol(tuples)){
     ## handles non-uniform arrays (cands depends on tuple)
     if (!uniform){
       cands <- do.call(expand.grid, mapply(":", rep(vstart,t), vend[tuples[,i]], SIMPLIFY = FALSE))
       ncands <- nrow(cands)
     }
     ticked <- rep(FALSE, nrow(cands))
      for (j in 1:N){
        nowrow <- D[j,tuples[,i]]
        if (flexval %in% nowrow) next   ## skips tuples that already hold an NA value
        now <- matrix(D[j,tuples[,i]], ncands, t, byrow=TRUE) ## t-vector of combinations to be covered
        nowcov <- which(apply(cands==now, 1, all))  ## or anyDuplicated(rbind(nowcov, cands)) - 1 better?
        if (!ticked[nowcov]){
          hilf[j, tuples[,i]] <- FALSE
          ticked[nowcov] <- TRUE
        }
      }
   }
   aus <- D; aus[hilf] <- NA
   ## profile not calculated with hilf, because D might already have NA entries
   attr(aus, "flexible") <- list(value=NA,
                                 profile=flexprofile(aus))
   aus
}

#' @export
flexpos <- function(D, t, ...){
  ## creates logical matrix with flexible positions set to TRUE
  if (any(is.na(D))) message("Note: Some flexible positions are already fixed. This reduces the positions that can be found.")

  hilf <- matcheck(D, uniform=FALSE, flexible=NA)
  v <- hilf$v ## vector of length k;
  k <- hilf$k; N <- hilf$N
  uniform <- FALSE
  if (length(unique(v))==1){
    uniform <- TRUE
    v <- v[1]
  }

  start0 <- hilf$start0
  vstart <- as.numeric(!start0)   ## number
  vend <- v - as.numeric(start0)  ## number or vector

  flexval <- NA
  ## start from all positions being redundant
  flexposition <- matrix(TRUE, N, k)
  tuples <- nchoosek(k,t)
  for (i in 1:ncol(tuples)){
    ## find duplicated entries
    now <- D[,tuples[,i]]
    ## incomparables does not work (would have liked to declare NA incomparable)
    dups <- union(which(duplicated(now)),
                       which(duplicated(now, fromLast=TRUE)))
    needed <- setdiff(1:N, dups)
    ## make non-duplicates non-flexible
    flexposition[needed, tuples[,i]] <- FALSE
  }
  flexposition[which(is.na(D))] <- TRUE
  flexposition
}

#' @export
flexprofile <- function(D, ...){
  hilf <- matcheck(D, uniform=FALSE, flexible=NA)
  v <- hilf$v ## vector of length k;
  k <- hilf$k; N <- hilf$N
  uniform <- FALSE
  if (length(unique(v))==1){
    uniform <- TRUE
    v <- v[1]
  }
  colSums(apply(D, 2, is.na))
}
