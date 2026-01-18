#' Methods for class CAcoverage
#'
#' print method for class CAcoverage
#'
#' @param x an object of class \code{CAcoverage}, created by function \code{\link{coverage}} or \code{coverage_iter}
#' @param ... further arguments to print
#' @param verbose logical, if TRUE, prints the verbose content as well
#'
#' @exportS3Method base::print
print.CAcoverage <- function(x, ..., verbose=FALSE){
   print(unlist(x[1:4]), ...)
   if (verbose) {
     cat("\n")
     print(x[-(1:4)], ...)
   }
}

