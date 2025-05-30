#' Methods for class coverage
#'
#' print method for class coverage
#'
#' @param x an object of class coverage, created by function coverage
#' @param ... further arguments to print
#' @param verbose logical, if TRUE, prints the verbose content as well
#'
#' @exportS3Method base::print
print.coverage <- function(x, ..., verbose=FALSE){
   print(unlist(x[1:4]), ...)
   if (verbose) {
     cat("\n")
     print(x[-(1:4)], ...)
   }
}

