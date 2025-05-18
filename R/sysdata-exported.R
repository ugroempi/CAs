#' Support data of package CAs
#' @description
#'  This documentation file documents objects that ship in
#' sysdata.rda but are nevertheless exported, because they are of interest to
#' expert users. I have still not figured out how to get roxygen to export this.
#'
#' NULL
#' @rawRD
#' \code{CAEX_lineages} contains a list (labeled by strengths, currently only strength t=2)\cr
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
#' \code{PALEYcat}, \code{DWYERcat}, \code{CKRScat}, \code{WKScat}, \code{NISTcat} and \code{TJcat}
#' are data.frame objects (\code{PALEYcat}, \code{DWYERcat},
#' \code{CKRScat} and \code{WKScat}) or matrices with columns
#' \code{t}, \code{k}, \code{v} and \code{N}
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
#'  The Source for CAEX-related information is Torres-Jimenez et al. (2021),
#' the source for Meagher and Stevens related information is Meagher and Stevens (2005),
#' the source for LCDST cover starters is Lobb et al. (2012), the source for
#' \code{WKScat} is
#' Wagner et al. (2021), the source for \code{CKRScat} is Colbourn et al. (2010),
#' and the source for \code{PALEYcat} is Colbourn (2015).
#'
#' The original web sources for both the Colbourn catalogue and the Torres-Jimenez
#' catalogue are currently unavailable; the November 2024 status of the Colbourn catalogue
#' is temporarily provided at \url{https://github.com/ugroempi/CAs/blob/main/ColbournTables.md}.
