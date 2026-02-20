funcheck <- function(test, tuple){
  ## test a vector of length k (permutation of the values 0 to k-1)
  ## tuple a vector of length t (sub sequence consisting of
  ##      t distinct ones from those values
  erg <- sapply(tuple, function(obj) which(test==obj))
  ## the sequence is covered, if erg is sorted
  all(erg==sort(erg))
}

## Alternative version that returns the redundant elements directly
identify_redundant_elements <- function(A, t) {
  stopifnot(is.matrix(A))
  k <- ncol(A)
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

can_reshuffle_to_cover <- function(row, useful_mask, seq) {
  ## Check if row contains all values from seq
  if (!all(seq %in% row)) return(FALSE)
  # works as planned
  # print(seq)
  # print(useful_mask)
  # print(row)

  ## Find positions where seq values currently are
  positions_of_seq_values <- sapply(seq, function(obj) which(row == obj))

  ## Can reshuffle if all elements in useful positions are in the correct order
  ## If at most one position is TRUE in useful_mask
  ## reshuffling is possible regardless of the order of the seq elements

  if (sum(useful_mask[positions_of_seq_values]) <= 1){
    aus <- TRUE
    attr(aus, "seq_useful_pos") <- which(useful_mask[positions_of_seq_values])
    attr(aus, "positions_of_seq_values") <- positions_of_seq_values
    return(aus)
  } else {
    ## more than one useful element in seq
    aus <- all(positions_of_seq_values[which(useful_mask[positions_of_seq_values])]==sort(positions_of_seq_values[which(useful_mask[positions_of_seq_values])]))
    if (!aus) return(aus) else{
      attr(aus, "seq_useful_pos") <- which(useful_mask[positions_of_seq_values])
      attr(aus, "positions_of_seq_values") <- positions_of_seq_values
      return(aus)
    }
  }
}

reshuffle_to_cover <- function(row, seq, positions_of_seq_values, seq_useful_pos) {
  ## This is a constraint satisfaction problem
  ## We need to rearrange non-useful elements so that seq appears in order
  ## This function should only be entered for cases with all elements of seq
  t <- length(seq)
  seq_moveable_pos <- setdiff(1:t, seq_useful_pos)
  ## all t values are redundant
  if (length(seq_moveable_pos) == t){
    new_row <- c(seq, setdiff(row, seq))
  }
  if (length(seq_moveable_pos) == t - 1){
    ## now only t-1 values are redundant
    before <- integer(0); after <- integer(0)
    if (seq_useful_pos > 1) before <- 1:(seq_useful_pos - 1)
    if (seq_useful_pos < length(seq)) after <- (seq_useful_pos + 1):length(seq)
    new_row <- c(seq[before], row[setdiff(1:length(row), setdiff(positions_of_seq_values, positions_of_seq_values[seq_useful_pos]))], seq[after])
  }
  ## now do the case with more than one useful element in seq
  if (length(seq_moveable_pos) < t - 1){
    moveable_pos_in_row <- positions_of_seq_values[seq_moveable_pos]
    row_reduced <- row[-moveable_pos_in_row]  # keep fixed seq elements and non-seq elements

    # Insert each moveable element just before the first fixed element
    # that must come after it in seq
    for (i in seq_along(seq_moveable_pos)) {
      seq_idx <- seq_moveable_pos[i]
      value_to_insert <- seq[seq_idx]

      useful_after <- seq_useful_pos[seq_useful_pos > seq_idx]
      if (length(useful_after) > 0) {
        anchor_value <- seq[min(useful_after)]
        anchor_pos <- which(row_reduced == anchor_value)
        row_reduced <- append(row_reduced, value_to_insert, after = anchor_pos - 1)
      } else {
        row_reduced <- c(row_reduced, value_to_insert)
      }
    }
    new_row <- row_reduced
  }
  new_row
}

## Helper: Actually reshuffle a row to cover a sequence
# reshuffle_to_cover <- function(row, seq, positions_of_seq_values, seq_useful_pos) {
#   ## This is a constraint satisfaction problem
#   ## We need to rearrange non-useful elements so that seq appears in order
#   ## This function should only be entered for cases with all elements of seq
#
#   t <- length(seq)
#   seq_moveable_pos <- setdiff(1:t, seq_useful_pos)
#   ## all t values are redundant
#   if (length(seq_moveable_pos) == t){
#     new_row <- c(seq, setdiff(row, seq))
#   }
#   if (length(seq_moveable_pos) == t - 1){
#     ## now only t-1 values are redundant
#     before <- integer(0); after <- integer(0)
#     if (seq_useful_pos > 1) before <- 1:(seq_useful_pos - 1)
#     if (seq_useful_pos < length(seq)) after <- (seq_useful_pos + 1):length(seq)
#     new_row <- c(seq[before], row[setdiff(1:length(row), setdiff(positions_of_seq_values, positions_of_seq_values[seq_useful_pos]))], seq[after])
#   }
#   ## now do the case with more than one useful element in seq
#   if (length(seq_moveable_pos) < t - 1){
#     new_row <- row
#     # Extrahiere die beweglichen Werte aus row
#     moveable_pos_in_row <- positions_of_seq_values[seq_moveable_pos]
#     moveable_values <- row[moveable_pos_in_row]
#
#     # Entferne die beweglichen Elemente aus row
#     row_reduced <- row[-moveable_pos_in_row]
#
#     # Bestimme für jedes bewegliche Element, wo es eingefügt werden muss
#     for (i in seq_along(seq_moveable_pos)) {
#       seq_idx <- seq_moveable_pos[i]
#       value_to_insert <- seq[seq_idx]
#
#       # Finde unbewegliche Elemente in seq, die VOR diesem Element kommen
#       useful_before <- seq_useful_pos[seq_useful_pos < seq_idx]
#
#       if (length(useful_before) > 0) {
#         # Finde das letzte unbewegliche Element in row_reduced
#         last_useful_value <- seq[max(useful_before)]
#         insert_after_pos <- which(row_reduced == last_useful_value)
#
#         # Füge danach ein
#         row_reduced <- append(row_reduced, value_to_insert, after = insert_after_pos)
#       } else {
#         # Ganz vorne einfügen
#         row_reduced <- c(value_to_insert, row_reduced)
#       }
#     }
#
#     new_row <- row_reduced
#   }
#   new_row
# }
