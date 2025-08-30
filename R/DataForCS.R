globalVariables(c("LCDSTStarters", "LCDSTCombis",
                  "CMMSSYStarters", "CMMSSYCombis",
                  "MeagherStevensStarters", "MeagherStevensCombis",
                  "ColbournKeriStarters", "ColbournKeriCombis"))

#' Data for cover starter constructions
#'
#' Lists ...Combis and ...Starters support cover starter constructions.
#'
#' @docType data
#'
#' @format
#' \code{LCDSTStarters}, \code{CMMSSYStarters}, \code{MeagherStevensStarters} and
#' \code{ColbournKeriStarters} are lists of lists:
#'
#' The first hierarchy level is either the number of levels per column for,
#' i.e., \code{v}, or the strength \code{t}:\cr
#' \code{CMMSSYStarters} holds 5 lists for 5 to 9 levels,\cr
#' \code{MeagherStevensStarters} holds 16 lists for 3 to 18 levels,\cr
#' \code{LCDSTStarters} 9 lists for 3 to 11 levels; these two constructions always yield strength 2.\cr
#' \code{ColbournKeriStarters} holds a single list for strength 4,
#' as the starter vectors for the other strengths can be calculated and do not need to be stored;
#' this construction works for 2 levels only.
#'
#' The second hierarchy level relates to the number of columns (\code{k}):\cr
#' Each list element at the first level is
#' a list whose entries are the starter vectors and whose names are either the
#' permissible numbers of columns or, for \code{CMMSSYStarters}, the length of the starter
#' which is one less than the number of columns.
#'
#' \code{MeagherStevensCombis} is a matrix with columns \code{v}, \code{k} and \code{N} that
#' indicates combinations for which starter vectors for the construction function \code{\link{CS_MS}} are available from
#' \code{MeagherStevensCombis}.
#'
#' \code{CMMSSYCombis} is a data.frame that
#' indicates combinations for which starter vectors for function \code{\link{CS_CMMSSY}} are available;
#' it has columns \code{v}, \code{k}, \code{N}, \code{k1}, \code{code}, and \code{stage}.\cr
#' \code{k1} indicates the number of columns in the left-hand side partition of a PCA structure.\cr
#' \code{code} is the code for creating the array, \code{stage} indicates the number of CAs
#' combined by function \code{\link{productPCA}}.
#'
#' \code{LCDSTCombis} is a data.frame that
#' indicates combinations for which starter vectors for function \code{\link{CS_LCDST}} are available;
#' it has columns \code{v}, \code{k}, \code{f}, \code{NAs}, \code{N},
#' \code{N_withoutInternet} and \code{N_eCANbased}.\cr
#' \code{f} indicates how many values are fixed and has an impact on the length of the starter vector;\cr
#' \code{N} is the number of runs achieved in the package based on installed constructions
#' for the ingredient CA with \code{f} levels, \code{N_withoutInternet} is that
#' same number if use of the internet is not possible, and \code{N_eCANbased}
#' is the number when using the best CA referenced in the Colbourn tables
#' (which may not be implemented in the package).
#'
#' \code{ColbournKeriCombis} is a data.frame with columns \code{t}, \code{k}, \code{v}, \code{N},
#' \code{method}, \code{starter}, \code{code} and \code{nconst} that
#' indicates combinations for which starter vectors for function \code{\link{CS_CK}} are available;
#' \code{method} refers to four different methods that are introduced by Colbourn and Keri (2009)
#' (one of which is deriving from the next higher strength),
#' \code{code} (built on \code{starter}) is used by function \code{CS_CK} for constructing
#' the actual arrays.
#'
#' @source
#' References are Meagher and Stevens (2005),
#' Lobb, Colbourn, Danziger, Stevens and Torres-Jimenez (2012, LCDST), and Colbourn and Keri (2009).
#'
#' @section Details:
#' For more details, see also the documentation of functions
#' \code{\link{CS_MS}}, \code{\link{CS_LCDST}} and \code{\link{CS_CK}}.\cr
#' The creation of the objects is documented by R files in the \code{extdata}
#' folder of this package, which can be located using
#' \code{system.file("extdata", package="CAs")}.
#'
#' @references Meagher and Stevens (2005), Colbourn and Keri (2009),
#' Lobb, Colbourn, Danziger, Stevens and Torres-Jimenez (2012, LCDST)
#'
#' @examples
#' ## an impression of the content
#' head(PALEYcat)
#' tail(PALEYcat)
#' boxplot(N ~ t, PALEYcat, las=1)
#'
#' head(CYCLOTOMYcat)
#' tail(CYCLOTOMYcat)
#' xtabs(~ t + v, CYCLOTOMYcat)
#'
#' ## the following is not run because of run times
#' ## in CRAN checks
#' \dontrun{
#' fun <- function(t,k,v) eCAN(t,k,v)$CAN
#' fun <- Vectorize(fun)
#' diff <- PALEYcat[,"N"]-fun(PALEYcat[,"t"], PALEYcat[,"k"], PALEYcat[,"v"])
#' fivenum(diff)
#' ## -102.0    2.0   35.5  168.0  438.0
#' ## one entry is much better than the current Colbourn table entry
#' ## this entry (t=6, k=380) has been checked for strength 6
#' ##     using samples of 20 columns each, as choose(380,6)
#' ##     was considered too large for a brute force check
#' ## no non-coverage was found, which is of course not proof
#' ##     for full strength 6 coverage
#' ##
#' }
#'

#'@rdname DataForCS
"MeagherStevensCombis"

#'@rdname DataForCS
"MeagherStevensStarters"

#'@rdname DataForCS
"LCDSTCombis"

#'@rdname DataForCS
"LCDSTStarters"

#'@rdname DataForCS
"CMMSSYCombis"

#'@rdname DataForCS
"CMMSSYStarters"

#'@rdname DataForCS
"ColbournKeriCombis"

#'@rdname DataForCS
"ColbournKeriStarters"
