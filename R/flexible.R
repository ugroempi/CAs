#' Redundancy checks and profile
#'
#' Functions for manipulating flexible values: determination of flexible positions (flexpos),
#' determination of a flexible set (markflex), and calculation of the profile (flexprofile)
#'
#' @rdname flexible
#'
#' @aliases flexpos
#' @aliases uniquecount
#' @aliases markflex
#' @aliases flexprofile
#' @aliases postopNCK
#'
#' @usage flexpos(D, t, ...)
#' @usage uniquecount(D, t, ...)
#' @usage markflex(D, t, fixrows=0, verbose=0, ...)
#' @usage flexprofile(D, ...)
#' @usage postopNCK(D, t, fixrows=0, verbose=0, outerRetry = 50, outerMaxnochange = 10,
#'        innerRetry=10, innerMaxnochange=25, seed=NULL, ...)
#'
#' @param D an N x k CA of strength t with v levels, coded from 0 to v-1 or from 1 to v; flexible values, if any, must be denoted as \code{NA}
#' @param t integer-valued (>=2), the strength of \code{D}
#' @param fixrows integer from 0 to \code{N} for preventing the top \code{fixrows} rows from being moved
#' @param verbose integer-valued degree of verbosity; 1 requests that the latest stages of unsuccessful attempts are added as an attribute to the returned unchanged matrix in case of a failure
#' @param outerRetry integer-valued number of retries for escaping local optimum, if \code{outerMaxnochange} iterations were not successful
#' @param outerMaxnochange integer-valued number of iterations to try and find optimum positions within a candidate row or its automatically-determined replacements
#' @param innerRetry integer-valued number of reshuffles for escaping local optimum, if \code{innerMaxnochange} iterations were not successful
#' @param innerMaxnochange integer-valued number of iterations to try and find optimum positions within a candidate row or its automatically-determined replacements
#' @param seed \code{NULL}, or an integer-valued seed for the random process;\cr
#'        the seed becomes an attribute of the output object
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
#' Function \code{flexprofile} obtains the profile of column-wise flexible values for \code{D}.
#'
#' Function \code{postopNCK} likely needs some tuning for speed and efficiency,
#' but is nevertheless often successful.\cr
#' It attempts to remove rows by trying to make entire rows
#' contain flexible values. It implements a variant of the proposed algorithm by Nayeri et al. (2013).
#' If successful, it returns \code{D} after removal of as many as possible flexible
#' rows, as created by \code{markflex}. If that immediate and cheap approach does not
#' work, the function tries to iteratively improve the number of flexible values in the
#' last row, until the row has flexible entries only. It that is not successful in a given
#' number of iterations (\code{innerMaxnochange}),
#' an inner retry step switches to a different last row to make flexible
#' (at most \code{innerRetry} attempts).\cr
#' An outer loop makes \code{outerRetry} such
#' attempts and keeps the user informed via messages about the current status.
#' If the user interrupts via the keyboard (\code{<Esc>} or \code{<Ctrl>-<c>}),
#' the result of the latest successful \emph{outer} retry is returned (make sure to
#' not abort iterations before the next outer retry was started).\cr
#' Without such an interrupt, if no attempt is successful,
#' the function returns the unchanged input matrix,
#' adding as an attribute a list with the last states of each retry,
#' if requested by the \code{verbose} argument.

#' @returns Function \code{flexpos} returns a logical matrix with \code{TRUE} values indicating
#' flexible positions of \code{D}. Any position that is already \code{NA} in \code{D} is set to \code{TRUE}.\cr
#' Function \code{markflex} returns the ingoing CA, with values denoted as before,
#' with rows potentially reordered (for facilitating all-flexible rows),
#' and with flexible positions changed to NA. \bold{this may change}\cr
#' Function \code{postopNCK} returns an object of class CA with prior attributes preserved
#' and attributes \code{Call} and \code{seed} augmented by the \code{postopNCK} call details.
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
#'
#' ## flexpos marks all positions that are not involved in any unique pairs
#' flexpos(plan, 2)
#' ## these cannot be simultaneously made flexible
#'
#' ## markflex creates actual flexible values (=don't care values)
#' ## and moves flexible rows to the end
#' (planMod <- markflex(plan, 2))
#' ## there are two rows with all values flexible (can be omitted)
#'
#' coverage(planMod[1:12,], 2)     ## now optimal
#' eCAN(2,7,3)
#'
#' ## fixing the first two rows from being moved
#' ## allows only one row to be removed
#' ## fixing rows comes at a price
#' (planMod2 <- markflex(plan, 2, fixrows=2))
#'
#' ## if removal of rows is the goal,
#' ## postopNCK does the entire job automatically
#' ## and very fast for this trivial case
#' plan_postopNCK <- postopNCK(plan, 2, outerRetry=1)
#' dim(plan_postopNCK)
#' coverage(plan_postopNCK, 2)
#'
#' plan_postopNCKfixtworows <- postopNCK(plan, 2, fixrows=2, innerRetry=2)
#' dim(plan_postopNCKfixtworows)
#' coverage(plan_postopNCKfixtworows, 2)
#'
#' ## postopNCK is fast on trivial problems like the above,
#' ## but may take its time otherwise

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
  suppressMessages( flexPositions <- flexpos(D, t, verbose=verbose, ...))
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
  attr(aus, "rowOrder") <- rowOrder
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
  ## perhaps change to iteration later
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


#' @export
postopNCK <- function(D, t, fixrows=0, verbose=0, outerRetry = 50, outerMaxnochange = 10,
                      innerRetry=10, innerMaxnochange=25, seed=NULL, ...){
  Call <- sys.call()
  attrs <- attributes(D)
  attrs$dim <- NULL
  levs <- levels.no.NA(D)
  ## lower bound for the run size
  bound <- prod(sort(levs, decreasing = TRUE)[1:2])
  if (nrow(D)==bound){
    message("D is already optimal")
    return(D)
  }
  Di <- markflex(D, t, fixrows=fixrows, verbose=verbose, ...)
  if (is.null(seed)) seed <- sample(32000,1)
  ncur <- nold <- nrow(D)
  count_unchanged <- 0
  set.seed(seed)
  message("seed: ", seed)

  ## obtained the technique via ChatGPT, prompts
  # In R, can I interrupt a function with Ctrl-C and obtain the solution of the latest completed iteration?
  # I want something that is permissible as a function in an R package;
  #     R packages are not permitted to write to the user's file system or the global environment
  #     without being specifically asked to do so. Is it possible to obtain the latest result in an
  #     assignment after an interrupt?
  last_result <- Di

  tryCatch({
  for (i in 1:outerRetry){
    message("postopNCK outer: ", i)
    ## call with default skipcheck=TRUE
    Di <- postopNCK_one(Di, t, fixrows=fixrows, verbose=verbose,
                        innerMaxnochange=innerMaxnochange, innerRetry=innerRetry, ...)
    if (ncur == nrow(Di)) {
      count_unchanged <- count_unchanged + 1
      if (count_unchanged >= outerMaxnochange) break
      }else{
      ncur <- nrow(Di)
      count_unchanged <- 0
      }
    message("run size: ", ncur)
    last_result <- Di
    if (ncur==bound){
      message("optimum is reached")
      break
      }
  }
  if (!"ca" %in% class(Di)) class(Di) <- c("ca", class(Di))
  attributes(Di) <- c(attributes(Di), attrs)
  attr(Di, "seed") <- c(attr(Di, "seed"), seed)
  attr(Di, "Call") <- Call
  return(Di)
  }, interrupt = function(e) {
    message("Interrupted by user. Returning last completed result.")
    if (!"ca" %in% class(last_result))
      class(last_result) <- c("ca", class(last_result))
    attributes(last_result) <- c(attributes(last_result), attrs)
    attr(last_result, "seed") <- seed
    attr(last_result, "Call") <- Call
    attr(last_result, "interrupt") <- i
    return(last_result)
  })
}

postopNCK_one <- function(D, t, fixrows=0, verbose=0, innerMaxnochange=25, innerRetry=10, skipcheck=TRUE, ...){
  if (!skipcheck) aus <- markflex(D, t, fixrows=fixrows, verbose=verbose, ...) else aus <- D
  orig_row <- attr(aus, "rowOrder")
  k <- ncol(aus)
  N <- nrow(aus)
  vs <- levels.no.NA(aus)
  hilf <- rowSums(is.na(aus))
  toremove <- which(hilf==k)
  if (length(toremove)>0){
    if (verbose>0) message("postopNCK removed ", length(toremove), " rows,\n",
                           N - length(toremove), " rows returned")
    orig_row <- orig_row[-toremove]
    aus <- aus[-toremove,]
    attr(aus, "rowOrder") <- orig_row
  }else{
    auslist <- vector(mode="list", length=innerRetry)
    exhausted <- rep(FALSE, N) ## initialize
    fixed <- numeric(0)
    if (fixrows > 0) {
      fixed <- 1:fixrows
      exhausted[fixed] <- TRUE
    }
    flexrows <- numeric(0)
    for (r in 1:innerRetry){
      message("inner retry: ", r)
      count <- 0
      pick <- which.max(hilf[!exhausted])
      ## positions must be moved, because exausted rows are prepended
      pick <- pick + sum(exhausted[1:pick])
      lastmax <- max(hilf)
      while (count <= innerMaxnochange){
      #if (hilf[pick] <= hilf[N]) pick <- N
      ## do the mixing, regardless whether moved or not
      mix <- sample(setdiff(1:N,c(fixed, pick)))
      exhausted <- exhausted[c(fixed, mix, pick)]
      orig_row <- orig_row[c(fixed, mix, pick)]
      aus <- aus[c(fixed, mix, pick), ]
      pick <- N ## now stays on this N for the entire r

      ## mixing rows helps to avoid stalls
      ## track exhausted row indicators along with matrix shuffles

      ## flexible positions of last row
      flexlast <- which(is.na(aus[N,]))
      ## put last row value into the flexible positions, where possible
      for (j in setdiff(1:k, flexlast)){
        if (any(is.na(aus[,j])))
        aus[which(is.na(aus[,j])),j] <- aus[N,j]
      }
      ## put random value into the flexible positions,
      ## where last row is flexible
        # assuming start0=TRUE

      for (j in flexlast){
        naj <- which(is.na(aus[,j]))
        aus[,j][naj] <- sample(0:(vs[j]-1), length(naj), replace=TRUE)
      }

      aus <- markflex(aus, t, fixrows=ifelse(!any(exhausted), fixrows, max(which(exhausted))), ...)
      orig_row <- orig_row[attr(aus, "rowOrder")]
      exhausted <- exhausted[attr(aus, "rowOrder")]
         if (exhausted[N]) cat("returned to exhausted last row ", orig_row[N],"\n")
      hilf <- rowSums(is.na(aus))
      curmax <- hilf[N] # max(hilf)
      pick <- N #which.max(hilf)
      if (curmax <=  lastmax){
        count <- count+1
        if (count==innerMaxnochange) cat(count, " changes did not increase beyond ", curmax, " flexible values\n")
        } else {
        lastmax <- curmax
        count <- 0
        }
      toremove <- which(hilf==k)
      if (length(toremove)>0){
        if (verbose>0) message("postopNCK removed ", length(toremove), " rows,\n",
                               N - length(toremove), " rows returned")
        aus <- aus[-toremove,]
        break
      }
    } ## end of while loop
      if (length(toremove) > 0) {
        attr(aus, "rowOrder") <- orig_row[-toremove]
        break ## return the aus from above
      }
      ## not yet successful, retry with different candidate row
      auslist[[r]] <- aus
      exhausted[N] <- TRUE
      flexrows <- setdiff(which(rowSums(is.na(aus)) > 0), fixed)
      aus <- aus[c(fixed, flexrows, setdiff((fixrows + 1):N,flexrows)),]
      orig_row <- orig_row[c(fixed, flexrows, setdiff((fixrows + 1):N,flexrows))]
      exhausted <- exhausted[c(fixed, flexrows, setdiff((fixrows + 1):N,flexrows))]

      ## replace missing values with random entries from 0 to vs[j]-1
      for (j in 1:k){
        aus[is.na(aus[,j]),j] <- sample(0:(vs[j]-1), sum(is.na(aus[,j])), replace=TRUE)
      }
      ## create new flexible set
      aus <- markflex(aus, t, fixrows=max(which(exhausted)))
      orig_row <- orig_row[attr(aus, "rowOrder")]
      exhausted <- exhausted[attr(aus, "rowOrder")]
      # print(orig_row[which(exhausted)])
      # print(orig_row[N])
      hilf <- rowSums(is.na(aus))
    }## end of loop over r
  }## end of else
  if (length(toremove)==0){
    if (verbose==1){
      attr(D, "postNCK_unsuccessful_attempts") <- auslist
    }
    return(D)
  }
  class(aus) <- c("ca", class(aus))
  if (any(is.na(aus)))
    attr(aus, "flexible") <- list(value=NA,
                                  profile=flexprofile(aus))
  aus
}


#' @export
uniquecount <- function(D, t, ...){
  ## creates logical matrix with flexible positions set to TRUE
  if (any(is.na(D))) {
    message("All flexible positions were fixed at 0.")
    D[is.na(D)] <- 0
  }

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
  uniquecount <- matrix(0, N, k)
  ## perhaps change to iteration later
  tuples <- nchoosek(k,t)
  for (i in 1:ncol(tuples)){
    ## find duplicated entries
    now <- D[,tuples[,i]]
    ## incomparables does not work (would have liked to declare NA incomparable)
    dups <- union(which(duplicated(now)),
                  which(duplicated(now, fromLast=TRUE)))
    needed <- setdiff(1:N, dups)
    ## make non-duplicates non-flexible
    uniquecount[needed, tuples[,i]] <- uniquecount[needed, tuples[,i]] + 1
  }
  uniquecount
}
