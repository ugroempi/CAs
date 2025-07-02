globalVariables("colbournBigFrame")

#' Colbourn tables for upper bounds on CAN / lower bounds on CAK
#'
#' data.frame that holds the November 2024 status of the Colbourn tables of upper bounds for covering array numbers
#'
#' @docType data
#'
#' @format
#' \code{data.frame} with 13641 rows in 5 columns:
#' \describe{
#'   \item{t}{the strength, from 2 to 6}
#'   \item{v}{the number of levels, from 2 to 25}
#'   \item{k}{the maximum number of columns}
#'   \item{N}{the smallest known run size for which a construction is *principally* available}
#'   \item{Source}{the table entry for the source information; sometimes (too) cryptic}
#' }
#'
#' @source  downloaded from \code{https://www.public.asu.edu/~ccolbou/src/tabby}
#' in November 2024; that source is no longer available, a slightly outdated version is available at
#' \url{https://web.archive.org/web/20231118015055/https://www.public.asu.edu/~ccolbou/src/tabby},
#' and a web version of the data frame is at \url{https://github.com/ugroempi/CAs/blob/main/ColbournTables.md}.\cr
#' The creation of the object is documented by the R file \code{colbournBigFrame_create.R} in the \code{extdata}
#' folder of this package, which can be located using
#' \code{system.file("extdata", package="CAs")}
#'
#' @section Details:
#' Until recently, the Colbourn (without year) tables were an authoritative source
#' of the state of the art run sizes of covering arrays and were kept reasonably
#' up-to-date by Charlie Colbourn. At present (June 2025), the tables are no longer
#' available to the public, and it is unclear, whether they are still updated.
#' The status of the tables in this package is November 2024. If the community finds
#' a way to maintain a database of this information, possibly in improved form,
#' this package needs to find a way to keep up with the changes. This
#' affects in particular the data.frame  \code{colbournBigFrame}, which is being used
#' as the basis behind decisions on which methods to implement.
#'
#' As the source information in the Colbourn tables is sometimes very terse,
#' it is not always intelligible for uninitiated humans,
#' including the package author.
#'
#' @references Colbourn (without year)
#'
#' @examples
#' head(colbournBigFrame)
#' ## Kleitman and Spencer () and Katona () provide the only construction
#' ## that is known to be globally optimal.
#' ## It is implemented in function \code{\link{KSK}}.
#'
#' tail(colbournBigFrame)
#' ## some entries have excessive N
#' ## this is also an example for the difficulty of deciphering the source entries:
#' ## the package author is not aware of the meanings of "N-CT" (what does the "N" stand for?),
#' ## "Arc", and "S", so that none of these constructions has been implemented.
#'
#' length(grep("projection", colbournBigFrame$Source))
#' length(setdiff(grep("projection", colbournBigFrame$Source),
#'                grep("postop", colbournBigFrame$Source)))
#' ## 137 of the 144 table entries that use a projection construction
#' ## were improved by postprocessing (which is not yet (?) available in this package)
#'

#'@rdname colbournBigFrame
"colbournBigFrame"

