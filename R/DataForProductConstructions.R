globalVariables(c("PCAcat", "DPcat"))

#' Data for productPCA and productCA constructions
#'
#' The data frames PCAcat and DPcat hold the information for strength 2 productPCA or productCA constructions
#' based on cover starter ingredients, Bose SCAs and additional good small CAs.
#'
#' @docType data
#'
#' @format
#' \code{PCAcat} is a data frame with columns \code{v}, \code{k} and \code{N}, \code{type},
#' \code{k1}, \code{code} and \code{stage}. \code{type} indicates which types of ingredients are
#' used, \code{k1} has the number of columns of the left-hand side parts in the resulting PCA,
#' \code{code} gives the code that is run for constructing the resulting CA, and \code{stage}
#' gives the number of ingredients multiplied for obtaining the result.
#'
#' \code{DPcat} is a data frame with columns \code{v}, \code{k} and \code{N}, \code{type},
#' \code{nconst}, \code{code} and \code{stage}. \code{type} indicates which types of ingredients are
#' used, \code{nconst} has the number of constant rows of the resulting PCA
#' (guaranteed, in a few cases slightly more would be possible),
#' \code{code} gives the code that is run for constructing the resulting CA, and \code{stage}
#' gives the number of ingredients involved in obtaining the result.
#'
#' @source
#' Colbourn and Torres-Jimenez (2013) and the references given below for the ingredients.
#'
#' @section Details:
#' The functions \code{\link{productPCA}} and \code{\link{productCA}} do the calculations,
#' based on ingredients created with the functions \code{\link{SCA_Bose}}, \code{\link{miscCA}},
#' \code{\link{CS_MS}}, \code{\link{CS_LCDST}} (\code{f=1} only) and \code{\link{CS_CMMSSY}}.\cr
#' The creation of the objects is documented by the R files \code{create_PCAcat.R} and
#' \code{create_DPcat.R} in the \code{extdata} folder of this package, which can be located using
#' \code{system.file("extdata", package="CAs")}. That folder also holds the candidate ingredients
#' used in the constructions in the file \code{PCA_cands.rda} and \code{DP_cands.rda}.
#'
#' @references Colbourn et al. (2006), Colbourn and Torres-Jimenez (2013), Bose (1938), Meagher and Stevens (2005), Colbourn and Keri (2009),
#' Lobb, Colbourn, Danziger, Stevens and Torres-Jimenez (2012, LCDST)
#'
#' @examples
#' ## an impression of the content
#' head(PCAcat)
#' tail(PCAcat)
#' boxplot(N ~ v, PCAcat, las=1, horizontal=TRUE)
#'
#' head(DPcat)
#' tail(DPcat)
#' boxplot(N ~ v, DPcat, las=1, horizontal=TRUE)
#'
#' ## the following is not run because of run times
#' ## in CRAN checks
#' \dontrun{
#' fun <- function(k,v) eCAN(2,k,v)$CAN
#' fun <- Vectorize(fun)
#' diff <- PCAcat[,"N"]-fun(PCAcat[,"k"], PCAcat[,"v"])
#' fivenum(diff)
#' ## 0    8   23  75  376
#' quot <- PCAcat[,"N"]/fun(PCAcat[,"k"], PCAcat[,"v"])
#' fivenum(quot)
#'
#' diff <- DPcat[,"N"]-fun(DPcat[,"k"], DPcat[,"v"])
#' fivenum(diff)
#' ## 0    8   23  75  376
#' quot <- DPcat[,"N"]/fun(DPcat[,"k"], DPcat[,"v"])
#' fivenum(quot)
#' }
#'

#'@rdname DataForProductConstructions
"PCAcat"

#'@rdname DataForProductConstructions
"DPcat"
