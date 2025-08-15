globalVariables(c("ca12.2.8", "ca26.4.15", "oa1728.12.6", "oa3375.15.6", "oa9261.21.6", "miscCAcat"))

#' Useful arrays for covering array constructions
#'
#' lists useful single arrays that are used themselves as covering arrays
#' or as ingredients for covering array constructions
#'
#' @docType data
#'
#' @format
#' All arrays listed here have their rows arranged such that constant rows are
#' at the top; if there are fewer than \code{v} constant rows, they are
#' arranged in PCA format (see \code{\link{productPCA}}).\cr
#' \code{ca26.4.15} is a strength 2 CA.\cr
#' \code{oa1728.12.6}, \code{3375.15.6} and \code{oa9261.21.6} are strength 3
#' orthogonal arrays (classes \code{ca} and \code{oa})
#' and matrices with levels starting at 0.
#' They might eventually reside in the R package \pkg{\link{DoE.base}}.\cr
#' \code{miscCAcat} is a data.frame that collects the properties of all the
#' CAs (including orthogonal arrays) that are documented in this file.
#' It has the columns \code{t}, \code{k}, \code{v},
#' \code{N}, \code{fns} (character), \code{nconst}, \code{PCAstatus} (numeric,
#' holds \code{k1}; 0 implies that \code{nconst=v}, i.e., k1=k)
#' and \code{hasNA} (logical, at present, FALSE for all entries).
#' .
#'
#' @source
#' \code{oa12.2.8} was created by the package author based on learning about
#' its existence (with two constant rows) from Colbourn and Torres-Jimenez (2010)
#' and Cohen, Colbourn and Ling (2003).\cr
#' \code{oa1728.12.6} was created by the package author using construction 2
#' of Ji and Yin (2010).\cr
#' The creation of the objects is documented by the R file \code{miscCAs.R}
#' in the \code{extdata} folder of this package, which can be located using
#' \code{system.file("extdata", package="CAs")}.
#'
#' @references Cohen, Colbourn and Ling (2003), Colbourn and Torres-Jimenez (2010), Ji and Yin (2010)
#'
#' @examples
#' attributes(ca12.2.8)
#' attributes(ca26.4.15)
#' attributes(oa1728.12.6)
#' attributes(oa3375.15.6)
#' attributes(oa9261.21.6)
#'

#'@rdname DataUsefulArrays
"ca12.2.8"

#'@rdname DataUsefulArrays
"ca26.4.15"

#'@rdname DataUsefulArrays
"oa1728.12.6"

#'@rdname DataUsefulArrays
"oa3375.15.6"

#'@rdname DataUsefulArrays
"oa9261.21.6"

#'@rdname DataUsefulArrays
"miscCAcat"
