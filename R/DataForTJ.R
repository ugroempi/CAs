globalVariables(c("TJcat", "TJ2level_CAs"))

#' Data for function tjCA
#'
#' The data.frame TJcat supports the identification of CAs from the
#' list TJ2level_CAs
#'
#' @docType data
#'
#' @format
#' \code{TJcat} is a data.frame with columns \code{t} (strength), \code{k} (number of CA columns),
#' \code{v} (number of levels), \code{N} (number of CA rows), \code{filesize} (the file size in bytes),
#' \code{nameInTJ2level_CAs}, \code{replaceable} and \code{code}.
#' Row names hold the names of the downloaded files.\cr
#' CAs corresponding to rows
#' with entries in \code{code} are created by executing the code of
#' column \code{code} (either with function \code{\link{CAEX}} for 3-level
#' or with one of several constructions for 2-level),\cr
#' CAs corresponding to rows with entries in \code{nameInTJ2level_CAs} by
#' picking the respective entry from \code{TJ2levelCAs}.\cr
#' CAs corresponding to rows with neither entry cannot be created at present
#' (many might be constructable by the power construction which will be
#' implemented shortly and will then be incorporated in the \code{code} column).
#'
#' \code{TJ2level_CAs} is a list of arrays (class \code{ca}) corresponding to the
#' rows of \code{TJcat} with entries in \code{nameInTJ2level_CAs}.
#'
#' @source
#' The source is the no longer functional website \code{https://www.tamps.cinvestav.mx/~oc/} by
#' Torres-Jimenez (without year). Web archive versions of that site are useless, because they
#' do not include the arrays.\cr
#' The creation of the objects is documented by R files in the \code{extdata}
#' folder of this package, which can be located using
#' \code{system.file("extdata", package="CAs")}, the actual raw arrays
#' can be found there as well (omitting large ones that can be constructed
#' at least approximately by a recursive construction).
#'
#' @section Details:
#' The data.frame \code{TJcat} holds information on 2-level CAs of strength 3
#' from Torres-Jimenez (2021),
#' covered by function \code{\link{CAEX}} of this package,\cr
#' and on 2-level CAs of strengths 3 to 6.
#' The list \code{TJ2level_CAs} holds actual arrays for the latter case.
#' These have been selected by size, i.e., too large ones have been omitted.
#' The unavailable 2-level arrays
#' (empty string in column \code{TJ2level_CAs}) can *in principle* be constructed
#' by the construction of Colbourn and Torres-Jimenez (2010).
#'  In various cases, post-processing steps by Torres-Jimenez, not implemented in the package,
#'  reduced the run size of the pure construction or added a column to the pure construction.
#'  Furthermore, as not all ingredients are available in this package (yet),
#' some sizes can not be achieved by the package. The size of the arrays
#' has nevertheless precluded the inclusion of the larger arrays into the package, in the interest
#' interest of limiting package size.
#'
#' @references Colbourn and Torres-Jimenez (2010), Torres-Jimenez (without year), Torres-Jimenez (2021)
#'
#' @examples
#' ## an impression of the content
#' head(TJcat)
#' fivenum(TJcat$N)
#' xtabs(~ v + t, TJcat)
#'
#' ## the following is not run because of run times
#' ## in CRAN checks
#' \dontrun{
#' fun <- function(t,k,v) eCAN(t,k,v)$CAN
#' fun <- Vectorize(fun)
#' diff <- TJcat[,"N"]-fun(TJcat[,"t"], TJcat[,"k"], TJcat[,"v"])
#' fivenum(diff)
#' ## -11 0 0 0 23
#' ratio <- TJcat[,"N"]/fun(TJcat[,"t"], TJcat[,"k"], TJcat[,"v"])
#' fivenum(ratio)
#' ## 0.9286 1.0000 1.0000 1.0000 1.0446
#' ## the arrays are pretty close to the best-known,
#' ##    occasionally even better
#' ##    (i.e., colbournBigFrame would have to be updated)
#' fun <- function(t,k,v) eCAN(t,k,v)$Source
#' fun <- Vectorize(fun)
#' }
#'

#'@rdname DataForTJ
"TJcat"

#'@rdname DataForTJ
"TJ2level_CAs"
