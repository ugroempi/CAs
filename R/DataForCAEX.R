globalVariables(c("CAEX_CAs", "CAEX_lineages"))

#' Data for function CAEX
#'
#' Two stored objects in support of function CAEX
#'
#' @docType data
#'
#' @format
#' \code{CAEX_CAs} is a list of stored CAs for the character string entries of
#' \code{CAEX_lineages}, with names corresponding to those character string entries.
#'
#' \code{CAEX_lineages} is a list of lists of lists of lists, where the first and second list
#' each have length one only and are named \code{2} for the strength and \code{3}
#' for the number of levels of all arrays, whereas the list \code{CAEX_lineages[['2']][['3']]}
#' has length 40 and names \code{11} to \code{50}, for the run sizes of the respective CAs.
#'
#' @source
#' The Source is Torres-Jimenez et al. (2021), and the arrays area also included in
#' \code{\link{TJcat}}.
#'
#' @section Details:
#' This section provides background information for interested users.
#' Array creation can be easily accomplished via the function \code{\link{CAEX}},
#' which works for strength 2 with 3-level columns.
#'
#' Torres-Jimenez et al. (2021) used product constructions and column extensions
#' to create a series of strength 2 CAs with 3-level columns from 11 to 50 runs.
#' Later CAs build on the best earlier CAs. The column extensions were obtained
#' with massive computer search effort. Hence, arrays obtained from those are stored
#' in the list \code{CAEX_CAs}, and their names are given in the object \code{CAEX_lineages}.
#' Arrays for larger numbers of runs (34, and 36 to 50) can be obtained from function
#' \code{\link{productCA}}, and the details are given in the respective entry of
#' \code{CAEX_lineage}.
#'
#' @references Torres-Jimenez (2021)
#'
#' @examples
#' lengths(CAEX_lineages[[1]][[1]])
#' CAEX_lineages[['2']][['3']][['12']]
#' CAEX_CAs[["ca12.2.3.7"]]
#' CAEX_lineages[['2']][['3']][['34']]
#'

#'@rdname DataForCAEX
"CAEX_lineages"

#'@rdname DataForCAEX
"CAEX_CAs"
