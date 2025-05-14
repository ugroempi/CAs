#' Internal data of package CAs
#'
#' This documentation file documents the internal objects that ship in
#' sysdata.rda. I have not yet found out how to make it part of the
#' official documentation, when using roxygen.
#'
#' @rdname internal-data
#' @docType data
#'
#' @aliases CAEX_CAs
#' @aliases CAEX_lineages
#' @aliases colbournBigFrame
#' @aliases colbournCatalogue
#' @aliases MeagherStevensCombis
#' @aliases MeagherStevensStarters
#' @aliases LCDSTCombis
#' @aliases LCDSTStarters
#' @aliases CKRScat
#' @aliases CKRS_CAs
#' @aliases WKScat
#' @aliases WKS_CAs
#' @aliases DWYERcat
#' @aliases NISTcat
#' @aliases TJcat
#'
#' @format \code{CAEX_lineages} contains a list (labeled by strengths, currently only strength t=2)\cr
#' of a list (labeled by number of levels, currently only v=3)\cr
#' of entries for each relevant run size N (labeled by this N, which ranges from 11 to 50 in the only implemented case);
#' list entries are character strings or lists of construction details.
#'
#' \code{CAEX_CAs} is a list of stored CAs for the character string entries of CAEX_lineages,
#' with names corresponding exactly to those character string entries.
#'
#' \code{WKS_CAs} is a list of stored strength 6 CAs for 2-level factors, constructed by Wagner, Kampel and Simos (2021).
#'
#' \code{CKRS_CAs} is a list of covering arrays provided on Keri's website, created according to Colbourn et al. (2010, CKRS)
#'
#' \code{colbournBigFrame} is a data.frame with columns \code{t}, \code{v}, code{k}, \code{N} and \code{Source}
#' of the content of the Colbourn tables,\cr
#' \code{colbournCatalogue} holds that same information as a list (for t=2,...,6) of lists (for v=2,...,25) of data.frame objects.
#'
#' \code{MeagherStevensCombis} is a matrix with columns \code{v}, \code{k} and \code{N} that
#' indicates combinations for which starter vectors for \code{\link{CS_MS}} are available,\cr
#' \code{MeagherStevensStarters} is a list of lists, with a list entry for each \code{v}
#' that holds a list of the starters for the relevant values of \code{k}.
#'
#' \code{LCDSTCombis} is a data.frame with columns \code{v}, \code{k}, \code{f}, \code{NAs} and \code{N} that
#' indicates combinations for which starter vectors for \code{\link{CS_LCDST}} are available,\cr
#' \code{LCDSTStarters} is a list of lists, with a list entry for each \code{v}
#' that holds a list of the starters for the relevant values of \code{k}.
#'
#' \code{DWYERcat}, \code{CKRScat}, \code{WKScat}, \code{NISTcat} and \code{TJcat} are data.frame objects (\code{DWYERcat},
#' \code{CKRScat} and \code{WKScat}) or matrices with columns \code{t}, \code{k}, \code{v} and \code{N}
#' (in different order, may change, always address via name!),
#' \code{DWYERcat} also with a \code{Source} column and \code{CKRScat} with a column \code{fn} for the file name
#' (that also indicates by containing x that the CA was constructed by cross-summing).\cr
#' They indicate the information on CAs that are available in the DWYER (2024) repository on Github,\cr
#' in this package downloaded from the Keri website accompanying the Colbourn et al. (2010) paper,\cr
#' in this package downloaded from the MATRIS website accompanying the Wagner et al. (2021) paper,\cr
#' in the NIST (2008) library of CAs
#' (not in the package, but publicly downloadable) or the Torres-Jimenez catalogue as of Feb 6, 2025
#' (currently non-public; planned to grow and make public, no timeline known).
#' The CAs in \code{CAEX_CAs} are from that catalogue.
#'
#' @source The Source for CAEX-related information is Torres-Jimenez et al. (2021),
#' the source for Meagher and Stevens related information is Meagher and Stevens (2005),
#' the source for LCDST cover starters is Lobb et al. (2012), the source for `WKScat` is
#' Wagner et al. (2021), and the source for `CKRScat` is Colbourn et al. (2010).
#'
#' The original web sources for both the Colbourn catalogue and the Torres-Jimenez
#' catalogue are currently unavailable; the November 2024 status of the Colbourn catalogue
#' is temporarily provided at \url{https://github.com/ugroempi/CAs/blob/main/ColbournTables.md}.
#'
