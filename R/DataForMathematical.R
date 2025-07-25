globalVariables(c("CYCLOTOMYcat", "PALEYcat"))

#' Data for function mathematical constructions
#'
#' The data.frames CYCLOTOMYcat and PALEYcat support the identification of valid
#' scenarii for use of the respective constructions, and the actual implementation.
#'
#' @docType data
#'
#' @format
#' \code{CYCLOTOMYcat} is a data.frame in 213 rows
#' with columns \code{t} (strength),
#' \code{k} (number of CA columns), \code{v} (number of levels), code{N} (number of CA rows),
#' \code{q} (prime or prime power behind the construcion), \code{type} (the construction type),
#' and the code for the workhorse function that the interface function \code{\link{cyclotomyCA}}
#' will use.
#'
#' \code{PALEYcat} is a data.frame in 178 rows with the same columns as
#' \code{CYCLOTOMYcat}, except for the column \code{type}, which is not needed;
#' here. \code{v} is always 2. The code is for the workhorse function
#' \code{paley_mat}, if necessary with rows added or removed, will be used by
#' the interface function \code{\link{paleyCA}}.
#'
#' @source
#' The reference for the cyclotomy construction is Colbourn (2010).
#' Some parameters in the paper were found to be wrong. Thus,
#' the parameters in \code{CYCLOTOMYcat} are from the Colbourn tables,
#' where hope is that they are correct.\cr
#' One instance from the Colbourn tables was removed,
#' as its apparent underlying q is not a prime,
#' but 13*229 so that either a different construction was used
#' or the entry carries wrong information.\cr
#' The cases that were checked
#' were found correct. It was not possible with the
#' available computing resources to check the actual coverage for
#' large cases, and sources apart from the Colbourn tables were neither
#' found for many of the arrays.
#'
#' The reference for Paley type constructions is Colbourn (2015),
#' and the constructions were taken from there.
#' For strength 3, only constructions
#' without doubling are included, as CAs from doubling are not competitive.
#' Small corrections to Colbourn's Table 1 were needed for primes 7,
#' 13 (strength 3 not possible) and 17 (no row can be removed).
#'
#' The creation of the objects is documented by R files in the \code{extdata}
#' folder of this package, which can be located using
#' \code{system.file("extdata", package="CAs")}.
#'
#' @section Details:
#' For more details, see also the documentation of functions
#' \code{\link{paleyCA}} and \code{\link{cyclotomyCA}}.
#'
#' @references Colbourn (2010), Colbourn (2015)
#'
#' @examples
#' ## an impression of the content
#' head(PALEYcat)
#' tail(PALEYcat)
#' boxplot(N ~ t, PALEYcat, las=1)
#'
#' head(CYCLOTOMYcat)
#' tail(CYCLOTOMYcat)
#' xtabs(~ t + v, CYCLOTOMYcat)
#'
#' ## the following is not run because of run times
#' ## in CRAN checks
#' \dontrun{
#' fun <- function(t,k,v) eCAN(t,k,v)$CAN
#' fun <- Vectorize(fun)
#' diff <- PALEYcat[,"N"]-fun(PALEYcat[,"t"], PALEYcat[,"k"], PALEYcat[,"v"])
#' fivenum(diff)
#' ## -102.0    2.0   35.5  168.0  438.0
#' ## one entry is much better than the current Colbourn table entry
#' ## this entry (t=6, k=380) has been checked for strength 6
#' ##     using 100 samples of 20 columns each, as choose(380,6)
#' ##     was considered too large for a complete brute force check.
#' ## no non-coverage was found, which is of course not proof
#' ##     for full strength 6 coverage.
#' ##
#' }
#'

#'@rdname DataForMathematical
"PALEYcat"

#'@rdname DataForMathematical
"CYCLOTOMYcat"
