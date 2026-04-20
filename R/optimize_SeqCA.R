#' Simulated Annealing for improving an SeqCA
#'
#' for k columns and strength t for the third stage of an algorithm by Torres-Jimenez et al. (2022)
#'
#' @rdname optimize_SeqCA
#'
#' @aliases optimize_SeqCA
#' @aliases simulated_annealing_sca
#'
#' @usage optimize_SeqCA(A, t, sa_params = list(), target=factorial(t),
#'            skipsa=FALSE, skipreduce=FALSE, verbose = FALSE, ...)
#' @usage simulated_annealing_sca(A, t, maxAccepted = 100, maxAttempts = 500,
#'            initialTemperature = 100, finalTemperature = 0.1,
#'            temperatureReduction = 0.95, probreshuffle=0, verbose = FALSE, ...)
#'
#' @param A an SeqCA of strength \code{t}, or for \code{simulated_annealing_sca}, a partial SeqCA for which coverage is to be completed
#' @param t the strength for which the coverage of all permutations of length t must be guaranteed for all \code{combn(k,t)} t-element subsets of 1,...,k (k the number of columns of A)
#' @param maxAccepted positive integer, maximum accepted perturbations for each temperature; default 100
#' @param maxAttempts positive integer, maximum attempts for each temperature; default 500
#' @param initialTemperature positive number, default 100
#' @param finalTemperature positive number, less than initialTemperature, default 0.1
#' @param temperatureReduction positive number, less than 1; the closer to 1, the slower the temperature reduction, and the more temperatures are run; default 0.95, corresponding to 135 temperatures with the other defaults
#' @param probreshuffle probability of a proper reshuffle in case that is possible; value between 0 and 1; positive values do not work yet.
#' @param skipsa logical, default FALSE; if TRUE, simulated annealing is suppressed
#' @param skipreduce logical, default FALSE; if TRUE, the initial iterative reduction step is suppressed (useful if the SeqCA is the result of such a reduction)
#' @param verbose logical; if TRUE, printing progress information
#' @param target number of rows at which no further improvement is attempted; defaults to the number of permutations needed to cover all orders for a single t-tuple, which is likely unachievable for most relevant problem sizes
#' @param sa_params list of parameters whose default values are to be changed; all parameters not in the list keep their default states; the possible paramters are \code{maxAccepted}, \code{maxAttempts},
#'  \code{initialTemperature}, \code{finalTemperature}, \code{temperatureReduction}, \code{probreshuffle}
#' @param ... currently not used
#'
#' @returns in case of success, a strength t SeqCA with k columns is returned; all rows contain the integers 1 to k.
#'
#' @section Details:
#' The function \code{simulated_annealing_sca} implements the simulated annealing
#' approach of Torres-Jimenez et al. (2022).
#' It can be applied directly to a partial SeqCA, but its main purpose is being called by function
#' \code{optimize_SeqCA} for reducing the number of rows of an existing strength t SeqCA (see below).#'
#' The perturbation function was not described in detail by Torres-Jimenez et al. 2022; they only
#' stated that the correct order for a currently non-covered randomly selected sequence is forced
#' for the row for which this action creates the fewest missings.\cr
#' This is implemented in the internal perturbation function within function
#' \code{simulated_annealing_sca} by simply reordering the sequence elements
#' (which might of course break existing coverages of the row, but might also bring useful new coverages).\cr
#' \code{optimize_SeqCA} implements a workflow for applying
#' \code{\link{reduce_rows_iterative_complete}} and subsequently applying simulated
#' annealing repeatedly.
#'
#' Torres-Jimenez et al. (2022) describe a 3-stage approach, the third stage of which is
#' simulated annealing. Simulated annealing is invoked if postprocessing of an initial greedy SeqCA
#' is unable to find further rows that can be removed without deteriorating coverage.
#' The post-processing and simulated annealing are combined in function \code{optimze_SeqCA}.
#' While the simulated annealing is needed for
#' achieving the very good performance that the paper reports on,
#' it can be very slow.
#' Therefore, it can be suppressed by setting \code{skipsa=TRUE}, in which case the function
#' \code{optimize_SeqCA} is basically a wrapper for \code{\link{reduce_rows_iterative_complete}}.\cr
#' The paper claims that the method is most useful for larger problems;
#' but for these, the run time will often be prohibitive.
#'
#' If \code{optimize_SeqCA} is stuck in simulated annealing iterations,
#' having achieved some progress but neither improving nor finishing, the process can
#' be interrupted by pressing the <ESC>-key, in which case it will return the last successful outcome.
#'
#' @section Use of AI:
#' Claude 4 was heavily involved in the development of these functions,
#' by initially implementing the simulated annealing pseudo code from the paper;
#' it worked reasonably well for everything except implementing the perturbation.
#'
#' @section Warning:
#' The implementation of simulated annealing is very preliminary,
#' and the \code{probreshuffle} argument does not work yet for positive probabilities.
#'
#' @seealso [greedySeqCA_TJ()] for the greedy method by
#' Torres-Jimenez et al. (2022) and
#' [reduce_rows_iterative_complete()] for the
#' postoptimization stage method of that paper that is used as the first stage of function
#' \code{optimize_SeqCA}
#'
#' @author Ulrike Groemping
#'
#' @references Torres-Jimenez et al. (2022)
#'
#' @examples
#' ## a small SeqCA with originally 13 rows that is post-optimized quickly to 9
#' A_TJ <- greedySeqCA_TJ(6,3, postopt = FALSE, seed=2323)
#' dim(A_TJ)
#'
#' A_TJ_reduced <- optimize_SeqCA(A_TJ, 3, target=9, verbose=TRUE)
#' A_TJ_redwoSimAnneal <- optimize_SeqCA(A_TJ, 3, skipsa=TRUE, verbose=TRUE)
#' coverageSeqCA(A_TJ_reduced, 3)
#' ## one less row is possible, but requires a very slow simulated
#' ## annealing optimization
#' A_TJ_reduced8 <- optimize_SeqCA(A_TJ, 3, target=8, verbose=TRUE)
#' \dontrun{
#' ## without setting the target, one would have to wait until
#' ## the simulated annealing for reducing to 7 runs failed
#'   optimize_SeqCA(A_TJ_reduced, 3, verbose=TRUE, skipreduce=TRUE)
#' }
#'
#' \dontrun{
#' ## much slower
#' nrow(A <- greedySeqCA_Kuhn(6, 4, seed=75))  ## 52 before, 41 after reduce
#' aus <- optimize_SeqCA(A, 4, skipreduce=TRUE, verbose=TRUE)
#'
#' }

#' @export
simulated_annealing_sca <- function(A, t,
                                    maxAccepted = 100,
                                    maxAttempts = 500,
                                    initialTemperature = 100,
                                    finalTemperature = 0.1,
                                    temperatureReduction = 0.95,
                                    probreshuffle = 0,
                                    verbose = FALSE, ...) {
  ## Input: A = partial SeqCA with some missing subsequences
  ## Output: complete or improved partial SeqCA


  k <- ncol(A)

  ## Generate all t-subsequences for reference
  combs <- combn(k, t)
  perms <- combinat::permn(t)
  # perms <- arrangements::permutations(n = t)
  # perms <- lapply(1:factorial(t), function(obj) perms[obj,])
  all_subsequences <- do.call(rbind, unlist(lapply(1:ncol(combs), function(obj)
    lapply(perms, function(obj2) combs[,obj][obj2])), recursive = FALSE))

  ## Helper: Count missing subsequences
  missings_of <- function(sca) {
    covered <- rep(FALSE, nrow(all_subsequences))
    for (i in 1:nrow(sca)) {
      for (seq_idx in 1:nrow(all_subsequences)) {
        if (!covered[seq_idx]) {
          if (funcheck(sca[i, ], all_subsequences[seq_idx, ])) {
            covered[seq_idx] <- TRUE
          }
        }
      }
    }
    sum(!covered)
  }

  ## Helper: Get list of non-covered subsequences
  get_non_covered <- function(sca) {
    covered <- rep(FALSE, nrow(all_subsequences))
    for (i in 1:nrow(sca)) {
      for (seq_idx in 1:nrow(all_subsequences)) {
        if (!covered[seq_idx]) {
          if (funcheck(sca[i, ], all_subsequences[seq_idx, ])) {
            covered[seq_idx] <- TRUE
          }
        }
      }
    }
    which(!covered)
  }

  ## Helper: Perturbation - modify SeqCA to try to cover a missing sequence
  perturbation <- function(sca, seq) {
    ## ## Strategy: Pick a random row and try to rearrange it to cover seq
    ## ## or swap elements between rows
    ## Strategy: determine the number of missings for each row after forcing seq to correct order
    ## keep row with fewest missings, randomly in case of ties
    new_sca <- sca
    n_rows <- nrow(sca)
    perturbed_missings <- rep(NA, n_rows)
    perturb_method <- rep("naive", n_rows)
    row_idx <- NULL
    for (i in 1:n_rows){
      new_sca <- sca
      positions <- which(new_sca[i,] %in% seq) ## sorted positions
      positions_in_order <- sapply(seq, function(obj) which(new_sca[i,] == obj))
      ## will useful stay the same?
      # if (can_reshuffle_to_cover(new_sca[i,], useful[i,], seq) && runif(1) < probreshuffle){
      #   if (funcheck(new_sca[i,], seq)) perturb_method[i] <- "dontchange" else{
      #   print(new_sca[i,])
      #   print(useful[i,])
      #   print(seq)
      #   print(positions_in_order)
      #    new_sca[i,] <- reshuffle_to_cover(new_sca[i,], seq, positions_in_order,
      #                                      which(useful[i,positions_in_order]))
      #    perturb_method[i] <- "reshuffle_to_cover"
      #   }
      # }
      # else
        new_sca[i, positions] <- seq
        ## perturb_method[i] remains unchanged

      ## regardless of method
      curmiss <- perturbed_missings[i] <- missings_of(new_sca)
      if (curmiss==0) {
        row_idx <- i
        perturbed_missings[i] <- 0
        if (verbose) cat("Found solution with 0 missings\n")
        break
      }
    }

    ## Pick a random row among the ones with the fewest missings
    if (is.null(row_idx)) row_idx <- which(perturbed_missings==min(perturbed_missings, na.rm = TRUE))
    if (length(row_idx) > 1) row_idx <- sample(row_idx, 1)
    new_sca <- sca

      ## Rearrange this row to cover seq
      ## Find positions of seq values
#    if (perturb_method[row_idx] == "naive"){
      positions <- which(new_sca[row_idx,] %in% seq)
      ## bring the sequence to its required order, if it is not yet

      ## If they're not in order, force them into the correct order
      ##    This could be made more sophisticated by using reshuffle...
      new_sca[row_idx, positions] <- seq
      # } else{
      #   if ()
      #   positions_in_order <- sapply(seq, function(obj) which(new_sca[i,] == obj))
      #   new_sca[i,] <- reshuffle_to_cover(new_sca[i,], seq,
      #                     positions_in_order, which(useful[i,positions_in_order]))
      # }
    new_sca
  }

  ## Step 3: Initialize temperature
  c <- initialTemperature

  ## Step 4: Initialize current solution
  currentSeqCA <- A
  currentMissings <- missings_of(A)

  ## Step 5: Initialize best solution
  bestSeqCA <- A
  bestMissings <- currentMissings

  if (verbose) {
    cat("Initial missings:", currentMissings, "\n")
    # cat("Target: 0 missings\n\n")
  }

  iteration <- 0

  ## Step 6: Main loop
  while (bestMissings > 0 && c > finalTemperature) {
    iteration <- iteration + 1
    if (verbose) cat("=== Temperature:", round(c, 3), ", Attempting missing removal for: ", nrow(currentSeqCA), "runs ===\n")

    ## Step 7: Reset counters
    accepted <- 0
    attempts <- 0

    ## Step 8: Inner loop
    while (bestMissings > 0 &&
           accepted <= maxAccepted &&
           attempts <= maxAttempts) {
      ## Step 9: Increment attempts
      attempts <- attempts + 1
      if (verbose && attempts%%10 == 0)
        cat("--- Temperature", round(c, 3), ", Attempt ", attempts, "---\n")

      ## Step 10: Select random non-covered subsequence
      non_covered_indices <- get_non_covered(currentSeqCA)
      if (length(non_covered_indices) == 0) break

      seq_idx <- sample(non_covered_indices, 1)
      seq <- all_subsequences[seq_idx, ]

      ## Step 11: Perturbation
      useful <- identify_needed(currentSeqCA, t)
      X <- perturbation(currentSeqCA, seq)

      ## Step 12: Calculate missings
      missings <- missings_of(X)

      ## Step 13: Acceptance criterion
      delta <- currentMissings - missings
      accept <- FALSE

      if (delta > 0) {
        ## Improvement - always accept
        accept <- TRUE
      } else {
        ## Worse or equal - accept with probability
        prob <- exp(delta / c)
        if (runif(1) < prob) {
          accept <- TRUE
        }
      }

      if (accept) {
        ## Step 14-15: Accept the change
        accepted <- accepted + 1
        currentSeqCA <- X
        currentMissings <- missings

        ## Step 16-18: Update best if improved
        if (currentMissings < bestMissings) {
          bestSeqCA <- currentSeqCA
          bestMissings <- currentMissings

          if (verbose) {
            cat("  New best: missings =", bestMissings,
                "(accepted =", accepted, ", attempts =", attempts, ")\n")
          }
        }
      }
    }

    if (verbose) {
      cat("  End of temperature level: accepted =", accepted,
          ", attempts =", attempts, ", best missings =", bestMissings, "\n\n")
    }

    ## Step 21: Reduce temperature
    c <- c * temperatureReduction
  }

  ## Step 23: Return best solution
  if (verbose) {
    cat("\n=== Final Result ===\n")
    cat("Best missings:", bestMissings, "\n")
    cat("Iterations:", iteration, "\n")
  }

  list(
    sca = bestSeqCA,
    missings = bestMissings,
    success = (bestMissings == 0),
    iterations = iteration
  )
}

#' @export
optimize_SeqCA <- function(A, t,
                         sa_params = list(),
                         target = factorial(t),
                         skipsa = FALSE,
                         skipreduce = FALSE,
                         verbose = FALSE, ...) {
  ## Helper: Count missing subsequences
  missings_of <- function(sca) {
    covered <- rep(FALSE, nrow(all_subsequences))
    for (i in 1:nrow(sca)) {
      for (seq_idx in 1:nrow(all_subsequences)) {
        if (!covered[seq_idx]) {
          if (funcheck(sca[i, ], all_subsequences[seq_idx, ])) {
            covered[seq_idx] <- TRUE
          }
        }
      }
    }
    sum(!covered)
  }

  stopifnot(coverageSeqCA(A, t)==1)
  if (nrow(A)<=target){
    message("A already has the target row size")
    return(A)
  }
  ## Start with initial SeqCA
  k <- ncol(A)
  current <- A

  if (verbose) {
    cat("=== Initial SeqCA ===\n")
    cat("Rows:", nrow(current), "\n\n")
  }
  ## start with run size reduction by reduce_rows_iterative
  if (!skipreduce){
  if (verbose)
    cat("=== Initial greedy reduction without simulated annealing ===\n")
  current <- reduce_rows_iterative_complete(current, t, verbose=verbose)
  if (nrow(current) <= target || skipsa){
    if (nrow(current) <= target && !skipsa)
      message("Initial reduction yielded target run size without simulated annealing")
    class(current) <- c("sca", class(current))
    return(current)
  }
  }

  if (verbose)
    cat("=== Prepare simulated annealing ===\n")

  if (!skipsa){
  ## Set default SA parameters
  sa_defaults <- list(
    maxAccepted = 100,
    maxAttempts = 500,
    initialTemperature = 100,
    finalTemperature = 0.1,
    temperatureReduction = 0.95,
    probreshuffle = 0
  )
  sa_params <- modifyList(sa_defaults, sa_params)

  iteration <- 0

  combs <- combn(k, t)
  perms <- combinat::permn(t)
  all_subsequences <- do.call(rbind, unlist(lapply(1:ncol(combs), function(obj)
    lapply(perms, function(obj2) combs[,obj][obj2])), recursive = FALSE))


  ## Iteratively remove rows and anneal to completeness
  last_result <- current

  tryCatch({

  repeat {
    iteration <- iteration + 1

    if (verbose) {
      cat("\n=== Iteration", iteration, "===\n")
      cat("Current rows:", nrow(current), "\n")
    }

  #   ## Identify useful elements
  #   useful <- identify_needed(current, t)
  #
    ## Try all possible row deletions and find the one with fewest missings
    best_deletion <- NULL
    best_deletion_missings <- Inf
    best_deletion_sca <- NULL

    for (i in 1:nrow(current)) {
      ## current is updated during the process
      ## therefore, this row is rerun at the beginning of each iteration step
      test_sca <- current[-i, ,drop=FALSE]
      missings <- missings_of(test_sca)
  #     result <- reduce_rows(current, useful, t, verbose = FALSE)
  #
  #     ## problem: i seems to change meaning during the process
  #     ## it is thus not guaranteed that all row deletions are tried (?)
  #
  #     # Check this specific row
  #     test_result <- list(
  #       deleted_row = i,
  #       missing_subsequences = 0
  #     )
  #
  #     # Actually test deleting row i
  #     test_sca <- current[-i, , drop = FALSE]
  #
  #     # Count missings
  #     covered <- rep(FALSE, nrow(all_subsequences))
  #     for (row_idx in 1:nrow(test_sca)) {
  #       for (seq_idx in 1:nrow(all_subsequences)) {
  #         if (!covered[seq_idx]) {
  #           if (funcheck(test_sca[row_idx, ], all_subsequences[seq_idx, ])) {
  #             covered[seq_idx] <- TRUE
  #           }
  #         }
  #       }
  #     }
  #     missings <- sum(!covered)
  #
      if (missings < best_deletion_missings) {
        best_deletion <- i
        best_deletion_missings <- missings
        best_deletion_sca <- test_sca
      }
      if (missings==0) break
    }

  #   if (is.null(best_deletion)) {
  #     if (verbose) cat("No row can be deleted\n")
  #     break
  #   }

    # useful <- identify_needed(current, t)
    # if (verbose) cat("=== remove a row with the fewest possible missings ===\n")
    # hilf <- reduce_rows(current, useful, t, verbose=verbose)
    # best_deletion_sca <- hilf$sca
    # if (verbose) {
    #   cat("Best deletion: row", hilf$deleted_row, "with", hilf$missings, "missings\n")
    # }

    ## If no missings, accept deletion and continue
    if (missings == 0) {
      current <- best_deletion_sca
      if (verbose) cat("Deletion creates no missings - accepted\n")
      next ## of repeat loop
    }

    ## If there are missings, try simulated annealing
    if (verbose) cat("Running simulated annealing to fix missings...\n")

    sa_result <- do.call(simulated_annealing_sca,
                         c(list(A = best_deletion_sca, t = t, verbose = verbose),
                           sa_params))

    ## Only accept if annealing succeeded (no missings)
    if (sa_result$success) {
      current <- sa_result$sca
      if (verbose) cat("Simulated annealing succeeded - deletion accepted\n")
      if (nrow(current)<=target) {
        if (verbose) cat("target run size reached")
        break
      }
    } else {
      if (verbose) {
        cat("Simulated annealing failed -", sa_result$missings, "missings remain\n")
        cat("Cannot reduce further - stopping\n")
      }
      break
    }
  }
  }, interrupt = function(e) {
    message("Interrupted by user. Returning last completed result.")
    if (!"sca" %in% class(last_result))
      class(last_result) <- c("sca", class(last_result))
    return(last_result)
  })
  }

  if (verbose) {
    cat("\n=== Final Optimized SeqCA ===\n")
    cat("Rows:", nrow(current), "\n")
    cat("Reduction from", nrow(A), "to", nrow(current), "rows\n")
  }
  class(current) <- c("sca", class(current))
  current
}
