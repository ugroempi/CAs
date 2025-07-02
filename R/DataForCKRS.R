globalVariables(c("CKRScat", "CKRS_CAs"))

#' Data for function ckrsCA
#'
#' The data.frame CKRScat supports the identification of CAs from the
#' list CKRS_CAs
#'
#' @docType data
#'
#' @format
#' \code{CKRScat} is a data.frame in 53 rows with columns \code{t} (strength), \code{v} (number of levels),
#' \code{k} (number of CA columns), code{N} (number of CA rows)
#' and \code{fn} (name of entry in list \code{CKRS_CAs}.
#'
#' \code{CKRS_CAs} is a list of 53 covering arrays (class \code{ca}) corresponding to the
#' rows of \code{CKRScat}.
#'
#' @source
#' The reference is Colbourn, Keri, Rivas Soriano and Schlage-Puchta (2010, CKRS).
#' The source for the arrays is \url{http://old.sztaki.hu/~keri/arrays/CA_listings.zip}.\cr
#' The creation of the objects is documented by R files in the \code{extdata}
#' folder of this package, which can be located using
#' \code{system.file("extdata", package="CAs")}, the actual raw arrays
#' can be found there as well (they are small enough to be included in spite
#' of their availability on the web).
#'
#' @section Details:
#' The data.frame \code{CKRScat} holds information on arrays in 2:7 and
#' 10:12 levels, at strengths 2 to 6.
#' The list \code{CKRS_CAs} holds the actual arrays.
#'
#' @references CKRS (2010)
#'
#' @examples
#' ## an impression of the content
#' head(CKRScat)
#' tail(CKRScat)
#' fivenum(CKRScat$N)
#' xtabs(~ v + t, CKRScat)
#'
#' ## the following is not run because of run times
#' ## in CRAN checks
#' \dontrun{
#' fun <- function(t,k,v) eCAN(t,k,v)$CAN
#' fun <- Vectorize(fun)
#' diff <- CKRScat[,"N"]-fun(CKRScat[,"t"], CKRScat[,"k"], CKRScat[,"v"])
#' fivenum(diff)
#' ## 0 0 2 6 69
#' ratio <- CKRScat[,"N"]/fun(CKRScat[,"t"], CKRScat[,"k"], CKRScat[,"v"])
#' fivenum(ratio)
#' ## 1 1 1.0010 1.0349 1.25
#' ## at least half of the arrays are pretty close to the best-known
#' }
#'

#'@rdname DataForCKRS
"CKRScat"

#'@rdname DataForCKRS
"CKRS_CAs"
