#' Reducing the number of rows for an SCA
#'
#' by covering all t-orders that are missing after removing a certain row through
#' rearrangements of redundant elements in other rows
#'
#' @rdname reduce_rows_complete
#'
#' @aliases reduce_rows_complete
#' @aliases reduce_rows_iterative_complete
#'
#' @usage reduce_rows_complete(A, t, verbose=FALSE, ...)
#' @usage reduce_rows_iterative_complete(A, t, verbose=FALSE, ...)
#'
#' @param A a sequence covering array (SCA) with k columns and elements 1,...,k in each row
#' @param t the strength for which the coverage of all permutations of length t must be guaranteed for all \code{combn(k,t)} t-element subsets of 1,...,k
#' @param verbose logical; if \code{TRUE}, printed output is more detailed
#' @param ... currently not used
#'
#' @returns an SCA with attributes
#'
#' @section Details:
#' This can likely be made faster by first ordering the rows in decreasing order of the count_effective attribute of coverage.
#'
#' @section Use of AI:
#' Claude 4 was heavily involved in the development of these functions.
#'
#' @examples
#' # a naive greedy SCA
#' A <- greedySCA_naive(6, 4, seed=17)
#' dim(A)  ## 67 rows
#' coverageSCA(A, 4) ## perfect coverage
#' \dontrun{
#'   Areduced <- reduce_rows_iterative_complete(A, 4)
#'   dim(Areduced)  ## 42 runs, i.e. reduced by about 1/3
#'   coverageSCA(Areduced,4)   ## still perfect coverage
#' }
#'

#'
#' @export
reduce_rows_complete <- function(A, t, verbose = FALSE, ...) {
  ## Input: A = SCA(N; t, k) matrix with perfect coverage
  ##        useful = matrix of same dimensions indicating useful elements
  ## Output: A with one row removed (if possible without losing coverage)
  stopifnot(coverageSCA(A, t) == 1)
  ## coverageSCA also checks the other properties of A

  Norig <- N <- nrow(A)

  useful <- identify_needed(A, t)

  redundant_rows <- which(rowSums(useful)==0)

  if (length(redundant_rows)>0){
    A <- A[-redundant_rows,]
    useful <- useful[-redundant_rows,]
    N <- nrow(A)
  }

  if (verbose && N < Norig){
    cat("immediate reduction from ", Norig, " to ", N, " rows\n")
    cat("by removing rows\n")
    print(redundant_rows)
  }

  redundant_rows <- integer(0)
  coverpct <- coverageSCA(A, t)
  all_subsequences <- attr(coverpct, "sequences")
  count_effective <- attr(coverpct, "count_effective")
  which_first <- attr(coverpct, "which_first")

  ## make sure that there are no "old" attributes
  # attributes(A)[c("deleted_row", "missings", "success")] <- NULL


  ## attempt further reductions
  ## Generate all t-subsequences for reference
     ## now taken from coverage
  # combs <- combn(k, t)
  # perms <- combinat::permn(t)
  # all_subsequences <- do.call(rbind, unlist(lapply(1:ncol(combs), function(obj)
  #   lapply(perms, function(obj2) combs[,obj][obj2])), recursive = FALSE))

  ## Step 3: Initialize
  deleteRow <- -1
  bestMissings <- Inf
  bestSCA <- A

  ## Step 4-22: For each row
  #for (i in N:1) { appears to work faster and better when starting at 1
  for (i in 1:N){
    ## outer loop: will be interrupted with break after a single success
    ## N as the upper loop bound does not change after starting the loop,
    ## though N will be reduced
     if (i > N) break
    if (verbose) cat("Testing removal of row", i, "\n")
    success <- FALSE  ## not proven that the row can be removed

    scaB <- A  ## Working copy
    usefulB <- identify_needed(scaB, t) ## Working copy
    which_firstB <- which_first ## Working copy that will be updated for currently tried i
                                ## which_first itself will be updated at successful removal of a row

    ## Step 6: Find subsequences involving useful elements of row i
    useful_in_row_i <- which(usefulB[i, ])
    stopifnot(length(useful_in_row_i) > 0) ## should not happen, because of immediate removal

    ## row does have useful elements
    ## Find which subsequences are covered by useful elements in row i
    seqs_to_check <- which_first[[i]]
    # for (seq_idx in 1:nrow(all_subsequences)) {
    #   ## inner loop over all subsequences
    #   seq <- all_subsequences[seq_idx, ]
    #   if (funcheck(A[i, ], seq)) {
    #     positions_used <- sapply(seq, function(val) which(A[i, ] == val))
    #     if (all(useful[i, positions_used])) {
    #       ## the sequence is covered for the first time in row i
    #       seqs_to_check <- c(seqs_to_check, seq_idx)
    #     }
    #   }
    # }

    if (verbose) cat("  Need to relocate", length(seqs_to_check), "subsequences\n")

    ## Step 6-13: For each subsequence involving useful elements
    for (ele in 1:length(seqs_to_check)) {
      seq_idx <- seqs_to_check[ele]
      seq <- all_subsequences[seq_idx, ]

      ## Step 7: Find another row (j != i) that already covers this sequence
      j <- NULL   ## replacement row for this seq

      for (row_j in setdiff(1:N, i)) {
        ## inner loop for replacement rows
        ################## this part should never be reached
        if (funcheck(scaB[row_j, ], seq)) {
            if (row_j < i) print("This is surprising and should not have happened")
            ## is already covered without reshuffling
            j <- row_j
            break ## break out of the loop over the runs
          }
        }
        if (!is.null(j)){
          ## make relevant elements useful
          usefulB[j, which(scaB[j,] %in% seq)] <- TRUE
          which_firstB[[j]] <- c(which_firstB[[j]], seq_idx)
          next  ## move to next sequence
        }
        ################# this part should never be run

      ## Step 8: If no row currently covers it, try to reshuffle
      if (is.null(j)) {
        ## Look for a row that COULD cover it by reshuffling non-useful elements
        for (row_j in setdiff(1:N, i)) {
          j_hilf <- can_reshuffle_to_cover(scaB[row_j, ], usefulB[row_j, ], seq)
          if (j_hilf) {
            j <- row_j
            seq_useful_pos <- attr(j_hilf, "seq_useful_pos")
            positions_of_seq_values <- attr(j_hilf, "positions_of_seq_values")
            break   ## break out of the loop over the runs
          }
        }
      }
      if (is.null(j)){
        ## no row found for this sequence
        break  ## from the loop over sequences
        ## next i will be gone to later
      } else {
        ## Step 8-12: Do the reshuffling
        scaB[j, ] <- reshuffle_to_cover(scaB[j, ], seq, positions_of_seq_values, seq_useful_pos)

          ## reshuffling for seq_idx may have led to new all-redundant row(s)
          usefulB[-i,] <- identify_needed(scaB[-i,], t)
          which_firstB[[j]] <- c(which_firstB[[j]], seq_idx)
   # print(funcheck(scaB[j,], seq))
          }
        } ## end of loop over to-be-checked sequences
        ## if got here without next, successful
        if (!is.null(j)) {
          success <- TRUE
          break
          } else {
            usefulB <- useful
            which_firstB <- which_first
            next
          }
      } ## end of loop over rows to be removed

    ## Step 14-18: Update best solution if this is better
    if (success) {
     # print("hier1")
      deleteRow <- i
      bestSCA <- scaB[-i, , drop = FALSE]
      useful <- identify_needed(bestSCA, t)
      which_first <- which_firstB[-i]
      if (any(rowSums(useful)==0)){
        bestSCA <- bestSCA[!rowSums(useful)==0,]
        useful <- useful[!rowSums(useful)==0,]
        which_first <- which_first[!rowSums(useful)==0]
      }
      N <- nrow(bestSCA)
      if (verbose) cat("  Achieved reduction to ", N, "rows from ", Norig, "rows\n")
    }
  return(bestSCA)
} ## end of function


#' @export
reduce_rows_iterative_complete <- function(A, t, verbose = FALSE, ...) {
  stopifnot(coverageSCA(A, t)==1) ## also checks properties of A
  original_rows <- nrow(A)
  k <- ncol(A)
  current_A <- A
  iter <- 0

  repeat {
    #useful <- identify_needed(current_A, t)
    result <- reduce_rows_complete(current_A, t, verbose = verbose)

    if (nrow(result) == nrow(current_A)) {
      if (verbose) cat("\nNo more rows can be removed\n")
      break
    }else{
      current_A <- result
      iter <- iter + 1
    }
  }

    if (verbose) {
      cat("\n=== Iteration", iter, "===\n")
      cat("Removed row(s), now have", nrow(current_A), "rows\n\n")
    }

  attr(current_A, "N_original") <- original_rows
  attr(current_A, "N_final") <- nrow(current_A)
  attr(current_A, "reduction_rate") <- (original_rows - nrow(current_A))/original_rows
  if (!"sca" %in% class(current_A)) class(current_A) <- c("sca", class(current_A))
  current_A
}
