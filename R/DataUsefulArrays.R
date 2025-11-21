globalVariables(c("ca12.2.8", "ca26.4.15", "oa1728.12.6", "oa3375.15.6", "oa9261.21.6",
                  "CohenSA", "ca23.4.8", "ca30.5.7", "ca41.6.6", "CohenSA", "miscCAcat"))

#' Useful arrays for covering array constructions
#'
#' lists useful single arrays that are used themselves as covering arrays
#' or as ingredients for covering array constructions
#'
#' @docType data
#'
#' @format
#' All arrays from this package that are listed here
#' have their rows arranged such that constant rows are
#' at the top or in PCA format; if there are fewer than \code{v} constant rows, they are
#' usually arranged in PCA format (see \code{\link{productPCA}}).\cr
#' \code{miscCAcat} also contains references to orthogonal arrays from R package \pkg{\link{DoE.base}},
#' e.g. \code{DoE.base::L81.3.5-1}, which is an OA(81,4,5,3)
#' (strength 4; 1 is subtracted, because OAs in that package
#' are coded with levels starting at 1,
#' whereas CAs in this package have levels starting at 0).\cr
#' \code{ca26.4.15}, \code{ca23.4.8}, \code{ca30.5.7} and \code{ca41.6.6} are strength 2 CAs.\cr
#' There is also the list \code{CohenSA} that contains strength 2 CAs that were created by Cohen using
#' simulated annealing.
#'
#' \code{oa1728.12.6}, \code{3375.15.6} and \code{oa9261.21.6} are strength 3
#' orthogonal arrays (classes \code{ca} and \code{oa})
#' and matrices with levels starting at 0.
#' They might eventually reside in the R package \pkg{\link{DoE.base}}.
#'
#' \code{miscCAcat} is a data.frame that collects the properties of all the
#' CAs (including orthogonal arrays and the CAs from within \code{CohenSA})
#' that are documented in this file.
#' It has the columns \code{t}, \code{k}, \code{v},
#' \code{N}, \code{fns} (character), \code{nconst}, \code{PCAstatus} (numeric,
#' holds \code{k1}; \code{v} implies that \code{nconst=v}, i.e., k1=k)
#' and \code{hasNA} (logical, at present, FALSE for all entries).
#' .
#'
#' @source
#' \code{oa12.2.8} was created by the package author based on learning about
#' its existence (with two constant rows) from Colbourn and Torres-Jimenez (2010)
#' and Cohen, Colbourn and Ling (2003).\cr
#' The list \code{CohenSA} holds arrays created by Cohen via simulated annealing (personal communication).\cr
#' \code{oa1728.12.6} was created by the package author using construction 2
#' of Ji and Yin (2010).\cr
#' The creation of the objects is documented by the R file \code{miscCAs.R}
#' in the \code{extdata} folder of this package, which can be located using
#' \code{system.file("extdata", package="CAs")}.
#'
#' @references Cohen, Colbourn and Ling (2003) and Cohen (personal communication), Colbourn and Torres-Jimenez (2010), Ji and Yin (2010), Kokkala et al. (2018)
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
"ca23.4.8"

#'@rdname DataUsefulArrays
"ca30.5.7"

#'@rdname DataUsefulArrays
"ca41.6.6"

#'@rdname DataUsefulArrays
"oa1728.12.6"

#'@rdname DataUsefulArrays
"oa3375.15.6"

#'@rdname DataUsefulArrays
"oa9261.21.6"

#'@rdname DataUsefulArrays
"CohenSA"

#'@rdname DataUsefulArrays
"miscCAcat"
