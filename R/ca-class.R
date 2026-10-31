#' Class ca and its methods
#'
#' This documentation file documents the class ca (for CAs = Covering Arrays) and its method.
#'
#' @rdname ca-class
#'
#' @aliases print.ca
#'
#' @method print ca
#'
#' @param x the object to be printed
#' @param \dots further arguments to \code{print.default}
#'
#' @returns The print method invisibly returns the input objext \code{x}
#'
#' @section Details:
#' The print method for class \code{ca} suppresses printing of most attributes,
#' but indicates which attributes are available.
#'
#' The S3 class \code{ca} is for a matrix with attributes.
#' Common attributes are \code{origin}, \code{constr}, \code{t}, \code{Call},
#' \code{eCAN}, \code{date}, \code{CAs-version}, \code{flexible}, \code{PCAstatus},
#' or \code{ingredients}. The \code{ingredients} attribute is a list of ingredient
#' constructions, each described by its \code{constr}, \code{t}, \code{k}, and \code{v}.
#' The \code{constr} attribute identifies the selected construction by its \code{Ns} label.
#' Not all attributes are present for all \code{ca} objects.
#'
#' The matrix elements are consecutive non-negative integers starting with 0 or 1,
#' and possibly missing values that stand for flexible values, i.e.,
#' values that can be arbitrarily chosen without deteriorating the strength of the CA.
#'
#' @export
print.ca <- function(x, ...){
   xnam <- deparse(substitute(x))
    if (!"ca" %in% class(x))
      stop("this print method is for class ca only")
   ## dim and dimnames should be kept by setdiff, i.e. not in attrs,
   ##    because dimnames is usually irrelevant and dim is needed
   ##    and dim is needed for printing as a matrix
    attrs <- setdiff(names(attributes(x)), c("origin", "class", "t", "dim", "dimnames"))
    for (a in attrs) attr(x, a) <- NULL
    print.default(x, ...)
    if (length(attrs) > 0) {
      message("further attribute(s)", "(accessible with attr(",
          xnam, ", attrname)):")
      print(attrs)
    }
    invisible(x)
 }
