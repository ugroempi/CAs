#' Identifying needed and redundant elements of SCAs
#'
#' from a given SCA A and its strength t
#'
#' @rdname identify_needed
#'
#' @aliases identify_needed
#'
#' @usage identify_needed(A, t, ...)
#'
#' @param A a sequence covering array (SCA) with k columns and elements 1,...,k in each row
#' @param t the strength for which the coverage of all permutations of length t must be guaranteed for all \code{combn(k,t)} t-element subsets of 1,...,k
#' @param ... currently not used
#'
#' @returns a logical matrix of the same dimensions as \code{A} with \code{TRUE} for needed values and \code{FALSE} for values that can be moved.
#'
#' @section The meanings of TRUE and FALSE in the resulting matrix:
#' In an SCA of strength t,
#' it is guaranteed that all orders of all t-tuples of elements of 1,...,k
#' are covered; the ordered tuples are also counted as covered if they are
#' interleaved by other values. For example, the sequence 1->4->3->2->5
#' covers the 3-sequence 1->3->5. The run 14325 covers \code{choose(5,3)}=10
#' 3-sequences, and each of its elements is part of six of these.
#' If all six such sequences are already covered in earlier runs,
#' the element is marked as redundant (\code{FALSE}),
#' otherwise it is needed (\code{TRUE}); for example, if the sequences 1->4->3, 1->4->2, 1->4->5,
#' 4->3->2, 4->3->5, 4->2->5 were already covered by earlier runs, the position of the "4" becomes \code{FALSE}.
#' Thus, if the position is \code{TRUE}, these sequences might still be covered
#' by later runs, but are not guaranteed to.
#'
#' @section Use of AI:
#' Claude 4 was involved in the development of this function.
#'
#' @examples
#' ### create a trivial SCA of strength 2
#' A <- rbind(1:5, 5:1)
#' identify_needed(A,2)
#'
#' ## a naive greedy SCA
#' A <- greedySCA_naive(5, 3, seed=2323, postop=FALSE)
#' dim(A)
#' ## identify the necessary elements
#' (neededs <- identify_needed(A,3))
#' ## there are various FALSE-only rows
#' ## remove those
#' A <- A[!rowSums(neededs)==0,]
#' neededs <- neededs[!rowSums(neededs)==0,]
#' dim(A)
#' ## all sequences are still covered
#' (cove <- coverageSCA(A, 3))
#' ## the count_effective attribute of cove shows
#' ## that the last three rows are only needed for 1 sequence each
#' attr(cove, "count_effective")
#' ## correspondingly, two elements of each of those rows are redundant
#' neededs
#' A[11:13,]
#' ## row 11 is needed for the sequence 3->1->4,
#' ## row 12 for the sequence 5->4->3,
#' ## and row 13 for the sequence 5->1->4
#'

#' @export
identify_needed <- function(A, t, ...) {
  ## Input: A = SCA(N; t, k) - a matrix with N rows and k columns
  ## Output: A matrix indicating which elements are needed (TRUE) or redundant (FALSE)
  ##         with the same dimensions as A

  stopifnot(is.matrix(A))
  k <- ncol(A)
  stopifnot(all(A %in% 1:k))

  ## Step 3: Set all elements in A as non-useful (FALSE)
  useful <- matrix(FALSE, nrow = nrow(A), ncol = ncol(A))

  ## Generate all possible t-subsequences
  ## These are all combinations of t positions from 1:k
  combs <- combn(k, t)
  # perms <- arrangements::permutations(n = t)
  # perms <- lapply(1:factorial(t), function(obj) perms[obj,])
  perms <- combinat::permn(t)

  ## All t-tuples to cover
  all_subsequences <- do.call(rbind, unlist(lapply(1:ncol(combs), function(obj)
    lapply(perms, function(obj2) combs[,obj][obj2])), recursive = FALSE))

  ## Step 4-7: For each subsequence
  for (seq_idx in 1:nrow(all_subsequences)) {
    seq <- all_subsequences[seq_idx, ]

    ## Step 5: Find first row that covers this subsequence
    for (i in 1:nrow(A)) {
      if (funcheck(A[i, ], seq)) {
        ## Step 6: Set the elements of seq in row i as useful
        ## Need to mark the positions in A[i,] that correspond to the values in seq
        for (val in seq) {
          pos <- which(A[i, ] == val)
          useful[i, pos] <- TRUE
        }
        break  ## Only the first row that covers it
      }
    }
  }

  ## Step 8: Return matrix indicating useful/redundant elements
  useful
}
