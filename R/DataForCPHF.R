globalVariables(c("WCS_CPHFs", "CL_SCPHFs", "CPHFcat"))

#' Data for the CPHF-based (permutation vector) constructions
#'
#' List WCS_CPHFs and data frame CPHFcat that support the CPHF-based constructions
#'
#' @docType data
#'
#' @format
#' \code{WCS_CPHFs} is a list of lists:
#'
#' The first hierarchy level is the strength \code{t}:\cr
#' There are four lists for strengths 3 to 6,\cr
#' each of which holds lists for 13 implemented numbers of levels \code{v},
#' which in turn hold lists of the implemented numbers of columns \code{k}.
#'
#' \code{CPHFcat} is a data frame with columns \code{t}, \code{v}, \code{k},
#' \code{nCPHF}, \code{hasNA} and \code{N}. The construction is from
#' Wagner, Colbourn and Simos (2022), with arrays provided by Michael Wagner.
#'
#' @source
#' Wagner, Colbourn and Simos (2022)
#'
#' @section Details:
#' For more details, see also the documentation of function
#' \code{\link{cphfCA}}.\cr
#' The creation of the objects is documented in the R file
#' \code{WCS_CPHFs_from_read.R} in the \code{extdata}
#' folder of this package, which can be located using
#' \code{system.file("extdata", package="CAs")};
#' the implemented raw CPHFs are not available in the package.
#'
#' @references Wagner, Colbourn and Simos (2022)
#'
#' @examples
#' ## an impression of the content
#' head(CPHFcat)
#' tail(CPHFcat)
#' names(WCS_CPHFs[["6"]][["3"]])
#' WCS_CPHFs[["6"]][["3"]][["7"]]  ## only one row --> CA with 3^6 runs
#' dim(WCS_CPHFs[["6"]][["3"]][["36"]])  ## 10 rows --> CA with 10*(3^6-1) + 1 runs
#'

#'@rdname DataForCPHF
"WCS_CPHFs"

#'@rdname DataForCPHF
"CPHFcat"
