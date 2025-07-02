globalVariables("NISTcat")

#' Data for function nistCA
#'
#' The data.frame NISTcat supports the identification and download of CAs from the
#' NIST library of covering arrays.
#'
#' @docType data
#'
#' @format
#' \code{NISTcat} is a matrix with columns \code{t} (strength), \code{v} (number of levels),
#' \code{k} (number of CA columns) and code{N} (number of CA rows). Its rownames are the file
#' names in the NIST library of covering arrays.
#'
#' @source
#' The source is NIST (2008).\cr
#' The creation of the object is documented by the R file \code{NISTcat_create.R}
#' in the \code{extdata} folder of this package, which can be located using
#' \code{system.file("extdata", package="CAs")}; the many
#' raw arrays are not provided there, as they are available on the web.
#'
#' @section Details:
#' The library is very large and has not been included into the package.
#' Function \code{\link{nistCA}} downloads the requested array as a zip file
#' into the temp directory (using function \code{\link{tempdir}} to locate
#' that directory), unzips it there and subsequently reads the CA into R.
#'
#' The CAs in the NIST library of covering arrays were obtained via computer
#' search. For many of them, there are meanwhile better alternatives.
#'
#' @references NIST (2008)
#'
#' @examples
#' ## an impression of the content
#' head(NISTcat)
#' tail(NISTcat)
#' xtabs(~t+v, as.data.frame(NISTcat))
#'
#' ## the following is not run because of run times
#' ## in CRAN checks
#' \dontrun{
#' fun <- function(t,k,v) eCAN(t,k,v)$CAN
#' fun <- Vectorize(fun)
#' diff <- NISTcat[,"N"]-fun(NISTcat[,"t"], NISTcat[,"k"], NISTcat[,"v"])
#' fivenum(diff)
#' ## 0 13 37 161 50698
#' ratio <- NISTcat[,"N"]/fun(NISTcat[,"t"], NISTcat[,"k"], NISTcat[,"v"])
#' fivenum(ratio)
#' ## 1 1.286 1.356 1.402 2.315
#' }
#'

#'@rdname DataForNIST
"NISTcat"
