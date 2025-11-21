globalVariables(c("SMC_SCPHFs", "CL_SCPHFs", "WC_SCPHFs", "SCPHFcat"))

#' Data for the SCPHF-based (permutation vector) constructions
#'
#' Lists SMC_SCPHFs, CL_SCPHFs and WC_SCPHFs and data frame SCPHFcat that support the SCPHF-based constructions
#'
#' @docType data
#'
#' @format
#' \code{SMC_SCPHFs}, \code{CL_SCPHFs} and \code{WC_SCPHFs}
#' are lists of lists:
#'
#' The first hierarchy level is the strength \code{t}:\cr
#' All lists hold at least elements for strengths 3 and 4, \code{WC_SCPHFs} also a few elements for strengths 5 and 6.\cr
#' Each of these holds lists for the implemented numbers of levels \code{v},
#' which in turn hold entries for the implemented numbers of columns \code{k}.
#'
#' \code{SCPHFcat} is a data frame with columns \code{t}, \code{v}, \code{k}, \code{N}
#' and \code{type}. The latter indicates the source: for "2006", the source is
#' the Sherwood et al. (2006) paper and SCPHFs are in \code{SMC_SCPHFs},\cr
#' for "2009" the source is Walker and Colbourn (2009) with arrays from Walker's
#' GitHub repository and arrays are in \code{WC_SCPHFs},\cr
#' for "2018" the source is work by Erin Lanus and co-authors,
#' and SCPHFs are in \code{CL_SCPHFs} (arrays found in the Dwyer (2024) GitHub repository).
#'
#' @source
#' Reference Sherwood, Martirosyan, Colbourn (2006), CPHFs and constructions
#' in the paper; Colbourn and Lanus (2018) and/or Colbourn, Lanus and Sarkar (2018),
#' SCPHFs from the repository of Dwyer (2024)
#'
#' @section Details:
#' For more details, see also the documentation of function
#' \code{\link{scphfCA}}.\cr
#' The creation of the objects is documented in the R files \code{SherwoodMartirosyanColbournCPHF.R},
#' \code{Lanus_SCPHFs_from_DwyerRepo.R} and \code{SCPHF_dimchecks.R} in the \code{extdata}
#' folder of this package, which can be located using
#' \code{system.file("extdata", package="CAs")}; the implemented raw SCPHFs from the DWYER
#' repo are in the sub folder \code{SCPHFs} of that folder.
#'
#' @references Colbourn and Lanus (2018), Colbourn, Lanus and Sarkar (2018), Dwyer (2024), Sherwood, Martirosyan, Colbourn (2006), Walker and Colbourn (2009)
#'
#' @examples
#' ## an impression of the content
#' SCPHFcat
#' SMC_SCPHFs[["4"]][["3"]]
#'

#'@rdname DataForSCPHF
"SMC_SCPHFs"

#'@rdname DataForSCPHF
"SCPHFcat"

#'@rdname DataForSCPHF
"CL_SCPHFs"

#'@rdname DataForSCPHF
"WC_SCPHFs"
