## Alternative version that returns the redundant elements directly
identify_redundant_elements <- function(A, t) {
  useful <- identify_needed(A, t)

  ## Return list with useful and redundant positions
  list(
    useful = useful,
    redundant = !useful,
    n_useful = sum(useful),
    n_redundant = sum(!useful),
    redundancy_rate = sum(!useful) / length(useful)
  )
}

## Visualization function
print_redundancy <- function(A, t) {
  result <- identify_redundant_elements(A, t)

  cat("SCA dimensions:", nrow(A), "x", ncol(A), "\n")
  cat("Total elements:", nrow(A) * ncol(A), "\n")
  cat("Useful elements:", result$n_useful, "\n")
  cat("Redundant elements:", result$n_redundant, "\n")
  cat("Redundancy rate:", sprintf("%.2f%%", result$redundancy_rate * 100), "\n\n")

  cat("Useful elements matrix:\n")
  print(result$useful)

  invisible(result)
}

## Example usage:
# A <- greedySCA_TJ(k=5, t=3)
# result <- identify_redundant_elements(A, t=3)
# print_redundancy(A, t=3)

reduce_rows <- function(A, useful, t, verbose = FALSE) {
  ## Input: A = SCA(N; t, k) matrix, could be partial only
  ##        useful = matrix of same dimensions indicating useful elements
  ## Output: A with one row removed (if possible without losing coverage)

  N <- nrow(A)
  k <- ncol(A)

  ## Generate all t-subsequences for reference
  combs <- combn(k, t)
  perms <- combinat::permn(t)
  # perms <- arrangements::permutations(n = t)
  # perms <- lapply(1:factorial(t), function(obj) perms[obj,])
  all_subsequences <- do.call(rbind, unlist(lapply(1:ncol(combs), function(obj)
    lapply(perms, function(obj2) combs[,obj][obj2])), recursive = FALSE))

  ## Step 3: Initialize
  deleteRow <- -1
  bestMissings <- Inf
  bestSCA <- A

  ## Step 4-22: For each row
  for (i in 1:N) {
    if (verbose) cat("Testing removal of row", i, "\n")

    scaB <- A  ## Working copy
    iMissings <- 0

    ## Step 6: Find subsequences involving useful elements of row i
    useful_in_row_i <- which(useful[i, ])

    if (length(useful_in_row_i) == 0) {
      ## No useful elements in this row - can delete without issues
      if (verbose) cat("  Row", i, "has no useful elements - can delete\n")
      deleteRow <- i
      bestSCA <- A[-i, , drop = FALSE]
      bestMissings <- 0
      break
    }

    ## Find which subsequences are covered by useful elements in row i
    seqs_to_check <- c()
    for (seq_idx in 1:nrow(all_subsequences)) {
      seq <- all_subsequences[seq_idx, ]
      if (funcheck(A[i, ], seq)) {
        positions_used <- sapply(seq, function(val) which(A[i, ] == val))
        if (all(useful[i, positions_used])) {
          seqs_to_check <- c(seqs_to_check, seq_idx)
        }
      }
    }

    if (verbose) cat("  Need to relocate", length(seqs_to_check), "subsequences\n")

    ## Step 6-13: For each subsequence involving useful elements
    for (seq_idx in seqs_to_check) {
      seq <- all_subsequences[seq_idx, ]

      ## Step 7: Find another row (j != i) that covers this sequence
      j <- NULL
      for (row_j in 1:N) {
        if (row_j != i && funcheck(scaB[row_j, ], seq)) {
          j <- row_j
          break
        }
      }

      ## Step 8: If no row currently covers it, try to reshuffle
      if (is.null(j)) {
        ## Look for a row that COULD cover it by reshuffling non-useful elements
        for (row_j in 1:N) {
          if (row_j == i) next

          if (can_reshuffle_to_cover(scaB[row_j, ], useful[row_j, ], seq)) {
            j <- row_j
            break
          }
        }
      }

      ## Step 8-12: Check if we found an alternative row
      if (!is.null(j)) {
        ## Step 9: Reshuffle non-useful elements in row j to cover seq
        ## the function needs the row, the sequence, the positions of sequence values,
        positions_of_seq_values <- sapply(seq, function(obj) which(scaB[j,]==obj))
        seq_useful_pos <- which(useful[row_j, positions_of_seq_values])
        if (!funcheck(scaB[j,], seq))
           scaB[j, ] <- reshuffle_to_cover(scaB[j, ], seq, positions_of_seq_values, seq_useful_pos)
      } else {
        ## Step 11: No alternative row found - this is a missing coverage
        iMissings <- iMissings + 1
      }
    }

    if (verbose) cat("  Removing row", i, "would leave", iMissings, "missing\n")

    ## Step 14-18: Update best solution if this is better
    if (iMissings < bestMissings) {
      deleteRow <- i
      bestSCA <- scaB[-i, , drop = FALSE]
      bestMissings <- iMissings

      if (verbose) cat("  New best: row", deleteRow, "with", bestMissings, "missing\n")
    }

    ## Step 19-21: If perfect solution found, return immediately
    if (bestMissings == 0) {
      if (verbose) cat("Perfect solution found - row", deleteRow, "can be removed\n")
      return(list(
        sca = bestSCA,
        deleted_row = deleteRow,
        missings = bestMissings,
        success = TRUE
      ))
    }
  }

  ## Step 23: Return best solution found
  if (verbose) {
    if (bestMissings == Inf) {
      cat("No row can be removed\n")
    } else {
      cat("Best option: remove row", deleteRow, "with", bestMissings, "missing subsequences\n")
    }
  }

  list(
    sca = bestSCA,
    deleted_row = deleteRow,
    missings = bestMissings,
    success = (bestMissings == 0)
  )
}


## Iterative version with complete implementation
reduce_rows_iterative <- function(A, t, verbose = FALSE) {
  original_rows <- nrow(A)
  k <- ncol(A)
  current_A <- A
  total_removed <- 0

  repeat {
    useful <- identify_needed(current_A, t)
    result <- reduce_rows(current_A, useful, t, verbose = verbose)

    if (!result$success || result$deleted_row == -1) {
      if (verbose) cat("\nNo more rows can be removed\n")
      break
    }

    total_removed <- total_removed + 1
    current_A <- result$sca

    if (verbose) {
      cat("\n=== Iteration", total_removed, "===\n")
      cat("Removed row, now have", nrow(current_A), "rows\n\n")
    }
  }

  list(
    sca = current_A,
    original_rows = original_rows,
    final_rows = nrow(current_A),
    rows_removed = total_removed,
    reduction_rate = total_removed / original_rows
  )
}

# final <- reduce_rows_iterative(A, t=3, verbose=TRUE)
# cat("\nFinal SCA:", nrow(final$sca), "rows\n")
# cat("Reduced from", final$original_rows, "to", final$final_rows, "rows\n")
# cat("Reduction:", sprintf("%.1f%%", final$reduction_rate * 100), "\n")

simulated_annealing_sca <- function(A, k, t,
                                    maxAccepted = 100,
                                    maxAttempts = 1000,
                                    initialTemperature = 100,
                                    finalTemperature = 0.1,
                                    temperatureReduction = 0.95,
                                    verbose = FALSE) {
  ## Input: A = partial SCA with some missing subsequences
  ## Output: complete or improved partial SCA

  ## Generate all t-subsequences for reference
  combs <- combn(k, t)
  perms <- combinat::permn(t)
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

## Complete workflow: Iterative row deletion with immediate annealing
optimize_sca <- function(A, k, t,
                         sa_params = list(),
                         verbose = FALSE) {
  ## Start with initial SCA
  current <- A

  if (verbose) {
    cat("=== Initial SCA ===\n")
    cat("Rows:", nrow(current), "\n\n")
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
      result <- reduce_rows_complete(current, useful, k, t, verbose = FALSE)

      # Check this specific row
      test_result <- list(
        deleted_row = i,
        missing_subsequences = 0
      )

      # Actually test deleting row i
      test_sca <- current[-i, , drop = FALSE]

      # Count missings
      combs <- combn(k, t)
      perms <- combinat::permn(t)
      all_subsequences <- do.call(rbind, unlist(lapply(1:ncol(combs), function(obj)
        lapply(perms, function(obj2) combs[,obj][obj2])), recursive = FALSE))

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

  current
}

# ## Example usage:
# A <- greedySCA_TJ(k=5, t=3)
# cat("Original:", nrow(A), "rows\n")
#
# # Optimize with default parameters
# optimized <- optimize_SCA(A, t=3, verbose=TRUE)
# cat("Optimized:", nrow(optimized), "rows\n")
#
# # Or just simulated annealing on a partial SCA
# partial <- A[1:10, ]  # Take first 10 rows (will have missings)
# result <- simulated_annealing_sca(partial, k=5, t=3, verbose=TRUE)
