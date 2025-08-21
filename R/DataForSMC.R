globalVariables(c("SMC_CPHFs", "SMCcat"))

#' Data for the SMC (permutation vector) construction
#'
#' Lists ...CPHFs and ...dat that support the SMC construction
#'
#' @docType data
#'
#' @format
#' \code{SMC_CPHFs} is a list of lists:
#'
#' The first hierarchy level is the strength \code{t}:\cr
#' \code{SMC_CPHFs} holds two lists for strengths 3 and 4,\cr
#' each of which holds lists for the implemented numbers of levels \code{v},
#' which in turn hold lists of the implemented numbers of columns \code{k}.
#'
#'
#' @source
#' Reference Sherwood, Martirosyan, Colbourn (2006), CPHFs and constructions
#' in the paper, one reduced from 40 to 32 columns because coverage is otherwise violated;
#' one omitted, because no simple fix of non-coverage was found
#'
#' @section Details:
#' For more details, see also the documentation of function
#' \code{\link{SMC}}.\cr
#' The creation of the object is documented by an R file in the \code{extdata}
#' folder of this package, which can be located using
#' \code{system.file("extdata", package="CAs")}.
#'
#' @references Sherwood, Martirosyan, Colbourn (2006)
#'
#' @examples
#' ## an impression of the content
#' SMCcat
#' SMC_CPHFs[["4"]][["3"]]
#'

#'@rdname DataForSMC
"SMC_CPHFs"

#'@rdname DataForSMC
"SMCcat"
