#' Checking coverage for SCAs
#'
#' for a given SCA A and a strength t
#'
#' @rdname coverageSCA
#'
#' @aliases coverageSCA
#'
#' @usage coverageSCA(A, t, stop_at_first_failure=FALSE, verbose=FALSE, ...)
#'
#' @param A a sequence covering array (SCA) with k columns and elements 1,...,k in each row
#' @param t the strength for which the coverage of all permutations of length t must be guaranteed for all \code{combn(k,t)} t-element subsets of 1,...,k
#' @param stop_at_first_failure logical; if \code{TRUE}, the function returns FALSE instead of a number
#' @param verbose logical; if \code{TRUE}, the returned object has an attribute \code{uncovered} with the uncovered sequences
#' @param ... currently not used
#'
#' @returns a number between 0 and 1 with attributes; the number is the proportion of covered orders (may be made more detailed in the future).\cr
#' The attribute \code{count_effective} holds the counts for each row of A of the number of sequences covered here but not earlier, as proposed by Murray and Colbourn (2015).\cr
#' If the proportion is smaller than one, and if \code{verbose=TRUE}, the uncovered sequences are available in the attribute \code{uncovered} of the returned value.
#'
#' @section Usefulness:
#' It is particularly relevant to check whether the coverage is perfect, i.e., the outcome is 1. If this is the only interest, one can stop at the first failure.
#'
#' @examples
#' ### create a trivial SCA of strength 2
#' A <- rbind(1:5, 5:1)
#' coverageSCA(A, 2)
#' coverageSCA(A, 3) ## covers a third of the 3-factor orders
#'
#' ## a naive greedy SCA of strength 3
#' A <- greedySCA_naive(5, 3, seed=2323)
#' dim(A)
#' coverageSCA(A,3)
#'
#' A <- reduce_rows_iterative_complete(A, 3)
#'
#' ## also works from worse starting matrix
#' A <- greedySCA_naive(5, 3, seed=589)
#' dim(A)
#' A <- reduce_rows_iterative_complete(A, 3)
#' dim(A)   ## substantially reduced in very short time
#'
#' @export
coverageSCA <- function(A, t, stop_at_first_failure=FALSE, verbose=FALSE, ...) {
  ## Input: A = SCA(N; t, k) - a matrix with N rows and k columns
  ## Output: a coverage proportion,
  ##         or the logical FALSE

  stopifnot(is.matrix(A))
  k <- ncol(A)
  stopifnot(all(A %in% 1:k))

  ## Generate all possible t-subsequences
  ## These are all combinations of t positions from 1:k
  combs <- combn(k, t)
  # perms <- arrangements::permutations(n = t)
  # perms <- lapply(1:factorial(t), function(obj) perms[obj,])
  perms <- combinat::permn(t)

  ## All t-tuples to cover
  all_subsequences <- do.call(rbind, unlist(lapply(1:ncol(combs), function(obj)
    lapply(perms, function(obj2) combs[,obj][obj2])), recursive = FALSE))
  stopifnot(nrow(all_subsequences) == (nseq <- choose(k,t)*factorial(t)))
  checkbits <- rep(0, nseq)  ## 0 for not covered, 1 for covered
  which_first <- vector(mode="list", length=nrow(A))

  ## Step 4-7: For each subsequence
  for (seq_idx in 1:nseq) {
    seq <- all_subsequences[seq_idx, ]

    ## Step 5: Find first row that covers this subsequence
    for (i in 1:nrow(A)) {
      if (funcheck(A[i, ], seq)) {
        checkbits[seq_idx] <- 1
        which_first[[i]] <- c(which_first[[i]], seq_idx)
        break  ## Only the first row that covers it
      }
    }
    if (stop_at_first_failure)
      if (checkbits[seq_idx]==0) stop("sequence ", paste(all_subsequences[seq_idx,], collapse="->"), " not covered")
  }
  aus <- mean(checkbits)
  attr(aus, "count_effective") <- lengths(which_first)
  attr(aus, "which_first") <- which_first
  attr(aus, "sequences") <- all_subsequences
  if (verbose && !stop_at_first_failure){
    uncovered <- character(0)
    if (any(checkbits==0))
    uncovered <- apply(all_subsequences[which(checkbits==0),,drop=FALSE], 1,
          function(obj) paste(obj, collapse="->"))
    attr(aus, "uncovered") <- uncovered
  }
  class(aus) <- c("coverageSCA", class(aus))
  aus
}

#' @export
print.coverageSCA <- function(x, ...){
  attributes(x) <- NULL
  print(x)
}
