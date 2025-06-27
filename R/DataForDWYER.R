globalVariables("DWYERcat")

#' Data for function dwyerCA
#'
#' The data.frame DWYERcat supports the identification and download of CAs from the
#' DWYER GitHub repository via function dwyerCA.
#'
#' @docType data
#'
#' @format
#' \code{DWYERcat} is a data.frame with columns \code{t} (strength), \code{v} (number of levels),
#' \code{k} (number of CA columns), code{N} (number of CA rows) and \code{Source}
#' (a brief descriptive source entry by Aaron Dwyer).
#'
#' @source
#' The source is Dwyer (2024), who cites personal communications with Colbourn as the
#' source of many of the arrays; there is some overlap between this repository and
#' \code{\link{CKRScat}}.
#'
#' @section Details:
#' The repository is large and has not been included into the package.
#' Function \code{\link{dwyerCA}} downloads calculates the requested array
#' (in case of an OA) or downloads it on the fly.
#'
#' @references Dwyer (2024)
#'
#' @examples
#' head(DWYERcat) ## will be constructed with function SCA_Busht
#' tail(DWYERcat) ## will be downloaded on the fly, except the last one
#'

#'@rdname DataForDWYER
"DWYERcat"
