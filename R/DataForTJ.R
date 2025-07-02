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
#' \code{v} (number of levels), code{N} (number of CA rows), \code{filesize} (the file size in bytes)
#' and \code{nameInTJ2level_CAs}. Row names hold the names of the downloaded files.
#'
#' \code{TJ2level_CAs} is a list of arrays (class \code{ca}) corresponding to the
#' rows of \code{TJcat} with element names.
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
#' by the construction of Colbourn and Torres-Jimenez (2010), which has not yet been
#' fully implemented. In variuos cases, post-processing steps reduced the run size of the
#' pure construction.
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
