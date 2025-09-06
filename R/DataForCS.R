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
#' \code{CMMSSYStarters} holds 12 lists for 3 to 10, 12, 13, 15, and 17 levels,\cr
#' \code{MeagherStevensStarters} holds 16 lists for 3 to 18 levels,\cr
#' \code{LCDSTStarters} 21 lists for 3 to 23 levels; these three constructions always yield strength 2.\cr
#' \code{ColbournKeriStarters} holds a single list for strength 4,
#' as the starter vectors for the other strengths can be calculated and do not need to be stored;
#' this construction works for 2 levels only.
#'
#' The second hierarchy level relates to the number of columns (\code{k}):\cr
#' Each list element at the first level is
#' a list whose entries are the starter vectors and whose names are the length of the starter;
#' except for \code{CMMSSYStarters}, the length of the starter equals the
#' permissible numbers of array columns, while for \code{CMMSSYStarters}, the array has
#' one more column than the length of the starter.
#'
#' \code{MeagherStevensCombis} is a matrix with columns \code{v}, \code{k} and \code{N} that
#' indicates combinations for which starter vectors for the construction function \code{\link{CS_MS}} are available from
#' \code{MeagherStevensCombis}.
#'
#' \code{CMMSSYCombis} is a data.frame that
#' indicates combinations for which starter vectors for function \code{\link{CS_CMMSSY}} are available;
#' it has columns \code{v}, \code{k}, \code{N}, \code{k1}, \code{code}, and \code{stage}.\cr
#' \code{k1} indicates the number of columns in the left-hand side partition of a PCA structure.
#'
#' \code{LCDSTCombis} is a data.frame that
#' indicates combinations for which starter vectors for function \code{\link{CS_LCDST}} are available;
#' it has columns \code{v}, \code{k}, \code{f}, \code{NAs}, \code{N},
#' \code{ingredient} and \code{N_eCANbased}.\cr
#' \code{f} indicates how many values are fixed;\cr
#' \code{N} is the number of runs achieved in the package based on installed constructions
#' for the ingredient CA with \code{f} levels,\cr
#' \code{ingredient} is the ingredient used in the construction (interesting for f>=4),\cr
#' and \code{N_eCANbased} is the number of runs achieved when using the best CA referenced in the
#' Colbourn tables (which may not be implemented in the package).
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
#' head(ColbournKeriCombis)
#' tail(ColbournKeriCombis)
#' length(ColbournKeriStarters[[1]])
#'
#' head(MeagherStevensCombis)
#' tail(MeagherStevensCombis)
#' length(MeagherStevensStarters)
#'
#' head(CMMSSYCombis)
#' tail(CMMSSYCombis)
#' length(CMMSSYStarters)
#'
#' head(LCDSTCombis)
#' tail(LCDSTCombis)
#' length(LCDSTStarters)
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
