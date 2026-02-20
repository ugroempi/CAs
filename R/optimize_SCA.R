#' Simulated Annealing for improving an SCA
#'
#' for k columns and strength t for the third stage of an algorithm by Torres-Jimenez et al. (2022)
#'
#' @rdname optimize_SCA
#'
#' @aliases optimize_SCA
#' @aliases simulated_annealing_sca
#'
#' @usage optimize_SCA(A, t, sa_params = list(), target=factorial(t), sa=TRUE, verbose = FALSE, ...)
#' @usage simulated_annealing_sca(A, t, maxAccepted = 100, maxAttempts = 1000, initialTemperature = 100, finalTemperature = 0.1, temperatureReduction = 0.95, verbose = FALSE, ...)
#'
#' @param A an SCA of strength \code{t}, or for \code{simulated_annealing_sca}, a partial SCA for which coverage is to be completed
#' @param t the strength for which the coverage of all permutations of length t must be guaranteed for all \code{combn(k,t)} t-element subsets of 1,...,k (k the number of columns of A)
#' @param maxAccepted positive integer
#' @param maxAttempts positive integer
#' @param initialTemperature positive number
#' @param finalTemperature positive number, less than initialTemperature
#' @param temperatureReduction positive number, less than 1; the closer to 1, the slower the temperature reduction, and the more iterations are run
#' @param sa logical, default TRUE; if FALSE, simulated annealing is suppressed
#' @param verbose logical; if TRUE, printing progress information
#' @param target number of rows at which no further improvement is attempted; defaults to the number of permutations needed to cover all orders for a single t-tuple, which is likely unachievable for most relevant problem sizes
#' @param sa_params list of parameters whose default values are to be changed; all parameters not in the list keep their default states
#' @param ... currently not used
#'
#' @returns in case of success, a strength t SCA with k columns is returned; all rows contain the integers 1 to k.
#'
#' @section Details:
#' The internal function \code{simulated_annealing_sca} implements the simulated annealing
#' approach of Torres-Jimenez et al. (2022) ??? check whether the perturbations are really the ones proposed there.
#' It is not meant to be applied directly but for being called by function
#' \code{optimize_SCA} for reducing the number of rows of an existing strength t SCA.
#' \code{optimize_SCA} implements a workflow for applying
#' \code{\link{reduce_rows_iterative_complete}} and subsequently applying simulated
#' annealing repeatedly.
#'
#' Torres-Jimenez et al. (2022) describe a 3-stage approach, the third stage of which is
#' simulated annealing. Simulated annealing is invoked if postprocessing of an initial greedy SCA
#' is unable to find further rows that can be removed without deteriorating coverage.
#' The post-processing and simulated annealing are combined in function \code{optimze_SCA}.
#' While the simulated annealing is needed for
#' achieving the very good performance that the paper reports on,
#' it is often excessively slow at least in the implementation of this package.
#' Therefore, it can be suppressed by setting \code{sa=FALSE}.\cr
#' The paper claims that the method is most useful for larger problems,
#' but for these, the run time will likely be prohibitive.
#'
#' @section Use of AI:
#' Claude 4 was heavily involved in the development of these functions,
#' by initially implementing the pseudo code from the paper.
#'
#' @seealso [greedySCA_TJ()] for the greedy method by
#' Torres-Jimenez et al. (2022) and
#' [reduce_rows_iterative_complete()] for the
#' postoptimization stage method of that paper
#'
#' @author Ulrike Groemping
#'
#' @references Torres-Jimenez et al. (2022)
#'
#' @examples
#' ## a small SCA with 13 rows that is optimized quickly to 9
#' A_TJ <- greedySCA_TJ(6,3, postopt = TRUE, seed=2323)
#' dim(A_TJ)
#'
#' A_TJ_reduced <- optimize_SCA(A_TJ, 3, target=9, verbose=TRUE)
#' A_TJ_redwoSimAnneal <- optimize_SCA(A_TJ, 3, sa=FALSE, verbose=TRUE)
#' coverageSCA(A_TJ_reduced, 3)
#' ## one less row is possible, but requires a very slow simulated
#' ## annealing optimization
#' \dontrun{
#' A_TJ_reduced8 <- optimize_SCA(A_TJ_reduced, 3, target=8, verbose=TRUE)
#' ## without setting the target, one would have to wait until
#' ## the simulated annealing for reducing to 7 runs failed
#' }
#'
#' @export
simulated_annealing_sca <- function(A, t,
                                    maxAccepted = 100,
                                    maxAttempts = 1000,
                                    initialTemperature = 100,
                                    finalTemperature = 0.1,
                                    temperatureReduction = 0.95,
                                    verbose = FALSE, ...) {
  ## Input: A = partial SCA with some missing subsequences
  ## Output: complete or improved partial SCA

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

  ## Helper: Perturbation - modify SCA to try to cover a missing sequence
  perturbation <- function(sca, seq) {
    ## Strategy: Pick a random row and try to rearrange it to cover seq
    ## or swap elements between rows

    new_sca <- sca
    n_rows <- nrow(sca)

    ## Pick a random row
    row_idx <- sample(n_rows, 1)

    ## Check if this row contains all values from seq
    if (all(seq %in% new_sca[row_idx, ])) {
      ## Try to rearrange this row to cover seq
      ## Find positions of seq values
      positions <- sapply(seq, function(val) which(new_sca[row_idx, ] == val))

      ## If they're not in order, try to rearrange
      if (!all(positions == sort(positions))) {
        ## Simple perturbation: swap two random elements in this row
        pos1 <- sample(k, 1)
        pos2 <- sample(k, 1)
        temp <- new_sca[row_idx, pos1]
        new_sca[row_idx, pos1] <- new_sca[row_idx, pos2]
        new_sca[row_idx, pos2] <- temp
      }
    } else {
      ## Swap elements between two random rows to bring seq values together
      row2_idx <- sample(setdiff(1:n_rows, row_idx), 1)
      pos1 <- sample(k, 1)
      pos2 <- sample(k, 1)
      temp <- new_sca[row_idx, pos1]
      new_sca[row_idx, pos1] <- new_sca[row2_idx, pos2]
      new_sca[row2_idx, pos2] <- temp
    }

    new_sca
  }

  ## Step 3: Initialize temperature
  c <- initialTemperature

  ## Step 4: Initialize current solution
  currentSCA <- A
  currentMissings <- missings_of(A)

  ## Step 5: Initialize best solution
  bestSCA <- A
  bestMissings <- currentMissings

  if (verbose) {
    cat("Initial missings:", currentMissings, "\n")
    cat("Target: 0 missings\n\n")
  }

  iteration <- 0

  ## Step 6: Main loop
  while (bestMissings > 0 && c > finalTemperature) {
    iteration <- iteration + 1
    if (verbose) cat("=== Temperature:", round(c, 3), "===\n")

    ## Step 7: Reset counters
    accepted <- 0
    attempts <- 0

    ## Step 8: Inner loop
    while (bestMissings > 0 && accepted <= maxAccepted && attempts <= maxAttempts) {
      ## Step 9: Increment attempts
      attempts <- attempts + 1

      ## Step 10: Select random non-covered subsequence
      non_covered_indices <- get_non_covered(currentSCA)
      if (length(non_covered_indices) == 0) break

      seq_idx <- sample(non_covered_indices, 1)
      seq <- all_subsequences[seq_idx, ]

      ## Step 11: Perturbation
      X <- perturbation(currentSCA, seq)

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
        currentSCA <- X
        currentMissings <- missings

        ## Step 16-18: Update best if improved
        if (currentMissings < bestMissings) {
          bestSCA <- currentSCA
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
    sca = bestSCA,
    missings = bestMissings,
    success = (bestMissings == 0),
    iterations = iteration
  )
}

#' @export
optimize_SCA <- function(A, t,
                         sa_params = list(),
                         target = factorial(t),
                         sa = TRUE,
                         verbose = FALSE, ...) {
  stopifnot(coverageSCA(A, t)==1)
  if (nrow(A)<=target){
    message("A already has the target row size")
    return(A)
  }
  ## Start with initial SCA
  k <- ncol(A)
  current <- A

  if (verbose) {
    cat("=== Initial SCA ===\n")
    cat("Rows:", nrow(current), "\n\n")
  }
  ## start with run size reduction by reduce_rows_iterative
  if (verbose)
    cat("=== Initial greedy reduction without simulated annealing ===\n")
  current <- reduce_rows_iterative_complete(current, t, verbose=verbose)
  if (nrow(current) <= target || !sa){
    if (nrow(current) <= target && sa)
      message("Initial reduction yielded target run size without simulated annealing")
    class(current) <- c("sca", class(current))
    return(current)
  }

  ## Set default SA parameters
  sa_defaults <- list(
    maxAccepted = 100,
    maxAttempts = 1000,
    initialTemperature = 100,
    finalTemperature = 0.1,
    temperatureReduction = 0.95
  )
  sa_params <- modifyList(sa_defaults, sa_params)

  iteration <- 0

  combs <- combn(k, t)
  perms <- combinat::permn(t)
  all_subsequences <- do.call(rbind, unlist(lapply(1:ncol(combs), function(obj)
    lapply(perms, function(obj2) combs[,obj][obj2])), recursive = FALSE))

  ## Iteratively remove rows and anneal to completeness
  repeat {
    iteration <- iteration + 1

    if (verbose) {
      cat("\n=== Iteration", iteration, "===\n")
      cat("Current rows:", nrow(current), "\n")
    }

    ## Identify useful elements
    useful <- identify_needed(current, t)

    ## Try all possible row deletions and find the one with fewest missings
    best_deletion <- NULL
    best_deletion_missings <- Inf
    best_deletion_sca <- NULL

    for (i in 1:nrow(current)) {
      ## current is updated during the process
      ## therefore, this row is rerun at the beginning of each iteration step
      result <- reduce_rows(current, useful, t, verbose = FALSE)

      ## problem: i seems to change meaning during the process
      ## it is thus not guaranteed that all row deletions are tried (?)

      # Check this specific row
      test_result <- list(
        deleted_row = i,
        missing_subsequences = 0
      )

      # Actually test deleting row i
      test_sca <- current[-i, , drop = FALSE]

      # Count missings
      covered <- rep(FALSE, nrow(all_subsequences))
      for (row_idx in 1:nrow(test_sca)) {
        for (seq_idx in 1:nrow(all_subsequences)) {
          if (!covered[seq_idx]) {
            if (funcheck(test_sca[row_idx, ], all_subsequences[seq_idx, ])) {
              covered[seq_idx] <- TRUE
            }
          }
        }
      }
      missings <- sum(!covered)

      if (missings < best_deletion_missings) {
        best_deletion <- i
        best_deletion_missings <- missings
        best_deletion_sca <- test_sca
      }
    }

    if (is.null(best_deletion)) {
      if (verbose) cat("No row can be deleted\n")
      break
    }

    if (verbose) {
      cat("Best deletion: row", best_deletion, "with", best_deletion_missings, "missings\n")
    }

    ## If no missings, accept deletion and continue
    if (best_deletion_missings == 0) {
      current <- best_deletion_sca
      if (verbose) cat("Deletion creates no missings - accepted\n")
      next
    }

    ## If there are missings, try simulated annealing
    if (verbose) cat("Running simulated annealing to fix missings...\n")

    sa_result <- do.call(simulated_annealing_sca,
                         c(list(A = best_deletion_sca, k = k, t = t, verbose = verbose),
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

  if (verbose) {
    cat("\n=== Final Optimized SCA ===\n")
    cat("Rows:", nrow(current), "\n")
    cat("Reduction from", nrow(A), "to", nrow(current), "rows\n")
  }
  class(current) <- c("sca", class(current))
  current
}
