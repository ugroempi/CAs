globalVariables(c("TJcat", "TJ2level_CAs"))

#' Data for function tjCA
#'
#' The data.frame TJcat supports the identification of CAs from the
#' list TJ2level_CAs or the creation of alternative arrays
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
#' rows of \code{TJcat} with entries in \code{nameInTJ2level_CAs} (access via names).
#'
#' @source
#' The source is the no longer functional website \code{https://www.tamps.cinvestav.mx/~oc/} by
#' Torres-Jimenez (without year), from which the arrays were downloaded on Feb 6th 2025.
#' Web archive versions of that site are useless, because they
#' do not include the arrays.\cr
#' The creation of the objects is documented by R files in the \code{extdata}
#' folder for this package, which can be located in the directory \code{inst/extdata}
#' of the GitHub repo.
#'
#' @section Details:
#' The data.frame \code{TJcat} holds information on 2-level CAs of strength 3
#' from Torres-Jimenez (2021),
#' covered by function \code{\link{CAEX}} of this package,\cr
#' and on 2-level CAs of strengths 3 to 6.
#' The list \code{TJ2level_CAs} holds actual arrays for the latter case.
#' These have been selected by size, i.e., too large ones have been omitted.
#' Where the unavailable 2-level arrays
#' (empty string in column \code{TJ2level_CAs}) can be replaced by other
#' constructions from the package, these are named in the \code{replaceable} column
#' (often the power construction with \code{powerCA} will yield such an array);
#' a \code{code} column provides the code for the replacement construction.
#'  Some sizes can not (yet) be achieved by the package. The size of the arrays
#' has nevertheless precluded the inclusion of the larger arrays into the package,
#' in the interest
#' interest of limiting package size.\cr
#' Note that a few replacement constructions also yield slightly smaller arrays
#' than those from the Torres-Jimenez repository.
#' Generally, users can simply use the \code{\link{Ns}} function or the \code{\link{bestN}} function to find
#' the best \bold{implemented} run size of this package; there is no need to work directly with
#' \code{TJcat} or \code{TJ2level_CAs}.
#'
#' @references Colbourn and Torres-Jimenez (2010), Torres-Jimenez (without year, download of arrays
#'  on Feb 6 2025), Torres-Jimenez (2021)
#'
#' @examples
#' ## an impression of the content
#' head(TJcat)
#' fivenum(TJcat$N)
#' xtabs(~ v + t, TJcat)
#'
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
#' ## constructions behind non-stored 2-level arrays
#' nonstored2level <- TJcat[TJcat$v==2 & TJcat$nameInTJ2level_CAs=="" & TJcat$replaceable=="", ]
#' table(substr(fun(nonstored2level[,"t"], nonstored2level[,"k"], nonstored2level[,"v"]), 1, 8))
#' ## the remaining not-yet-replaceable arrays are mostly doable from Derive and Power CZ
#' ## constructions, which have not systematically (Derive) or not at all (Power CZ) been implemented
#' ## in this package so far.

#'@rdname DataForTJ
"TJcat"

#'@rdname DataForTJ
"TJ2level_CAs"
