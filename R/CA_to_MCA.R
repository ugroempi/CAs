#' Making mixed-level CA from uniform CA
#'
#' takes a mixed level CA, makes unused levels flexible, and removes redundant rows
#'
#' @rdname CA_to_MCA
#'
#' @aliases CA_to_MCA
#'
#' @usage CA_to_MCA(D, cs, tolevs, t=attr(D, "t"), outerRetry=10, seed=NULL, ...)
#'
#' @param D a strength \code{t} CA (it may have mixed levels)
#' @param cs integer vector of column numbers
#' @param tolevs integer vector of target numbers of levels for the columns indicated
#'        by \code{cs}; must not be larger than the number of levels of the related
#'        columns in \code{cs}, and must not be smaller than 2
#' @param t integer, the strength of \code{D} (smaller strength would be permitted,
#'        but would not make sense)
#' @param outerRetry the number of calls to function \code{postopNCK}
#' @param seed the seed for making \code{\link{postopNCK}} reproducible;\cr
#'        if \code{NULL}, it will be randomly determined
#' @param ... currently not used
#'
#' @section Details:
#' The function sets unused levels of the columns in \code{cs} to \code{NA}, i.e.,
#' makes them flexible. Subsequently, the function \code{\link{postopNCK}} is used for
#' making entire rows flexible and removing them. Repeated runs can remove more rows.
#' The \code{outerRetry} argument gives the number of calls to \code{\link{postopNCK}}.
#' if the bottom rows of the outcome still look promising, the function can be called
#' again for removing even more rows.
#'
#' Without specifying a seed, the result will differ from call to call. The seed, even
#' if not explicitly specified, is stored with the output, so that the result will
#' always be repeatable, if desired.
#'
#' @returns a matrix of class \code{ca} with attributes of the ingoing matrix (dimensions modified),
#'    and attributes \code{Call} and \code{seed} modified.
#'
#' @author Ulrike Groemping
#'
#' @examples
#' # a small example
#' ## six 3-level columns in 33 runs
#' D <- bestCA(3,6,3)
#' ## to three 3-level and three 2-level columns in 27 runs
#' Dmixed <- CA_to_MCA(D, cs=4:6, tolev=rep(2, 3), t=3)
#' dim(Dmixed)
#' ## (optimal)
#' coverage(Dmixed, 3)
#' ## works
#'
#' ## five 3-level columns in 33 runs
#' D <- bestCA(3,5,3)
#' ## to four 3-level columns and one 2-level column in 33 runs,
#' ## i.e., a lengthy search cannot reduce the number of runs;
#' ## the 27 run OA(27,3,4,3) cannot easily be expanded by a 2-level column
#' ## without increasing the number of levels
#' \dontrun{
#' Dmixed <- CA_to_MCA(D, cs=5, tolev=2, t=3) ## last three columns reduced to 2 levels,
#' }
#'
#' ## larger example
#' D <- bestCA(3, 12, 5)
#' dim(D)
#' ## 225 runs for 555555555555 levels
#' ## can be reduced to 177 runs for 555554433222 levels
#' ## not run because of run time;
#' ## more reduction might be possible, though the eighth interim run
#' ## of the second round already failed
#' \dontrun{
#' Dmixed <- CA_to_MCA(D, cs=6:12, tolevs=c(4,4,3,3,2,2,2), seed=28661)
#' dim(Dmixed)  ## 184 x 12
#' head(Dmixed)
#' tail(Dmixed)
#' ## still looks promising for further row reductions
#' ## no problem to continue from the result array that already
#' ##     has the correct number of levels
#' ## could also use postopNCK directly, but no need
#' Dmixed2 <- CA_to_MCA(Dmixed, cs=6:12, tolevs=c(4,4,3,3,2,2,2), seed=7199)
#' dim(Dmixed2)
#' tail(Dmixed2)
#' ## there may still be potential for more reduction,
#' ## but not too bad
#' }
#'

#' @export
CA_to_MCA <- function(D, cs, tolevs, t=attr(D, "t"), outerRetry=10, seed=NULL, ...){
  Call <- sys.call()
  if (is.null(seed)) seed <- sample(32000,1)
  if (is.null(t) || !(is.numeric(t)))
    stop("D does not have a valid attribute 't',\nplease specify the strength t")
  stopifnot(is.matrix(D))
  stopifnot(is.numeric(D))
  start0 <- min(D)==0
  ## assuming valid array starting at 0 or 1
  k <- ncol(D)
  N <- nrow(D)
  levs <- levels.no.NA(D)
  stopifnot(is.numeric(cs))
  stopifnot(is.numeric(tolevs))
  stopifnot(all((c(t,cs,tolevs) %% 1) == 0))
  stopifnot(length(cs)==length(tolevs))
  stopifnot(all(tolevs>1))
  stopifnot(all(tolevs<=levs[cs]))
  attrs_stored <- attributes(D)
  levs[cs] <- tolevs
  optN <- lower_bound_CA(t, k, levs)
  for (j in 1:length(cs)){
    D[which(D[,cs[j]]>tolevs[j] - as.numeric(start0)),cs[j]] <- sample(levs[cs[j]]-as.numeric(start0),1)
  }
  if (optN == N){
    message("D has already the optimum number of rows.")
    attrs <- attrs_stored
    attrs$Call <- c(attrs_stored$Call, Call)
    attributes(D) <- attrs
    return(D)
  }
  Dnow <- D
  set.seed(seed)
  for (r in 1:outerRetry){
      aus <- postopNCK(Dnow, t)
      Dnow <- aus
      if (nrow(aus)==optN) break
  }
  if (nrow(aus) == N){
    message("The number of rows of D could not be reduced.")
    attrs <- attrs_stored
    attrs$Call <- c(attrs_stored$Call, Call)
    attributes(D) <- attrs
    return(D)
  }
  attr(aus, "rowOrder") <- NULL
  attrs <- attrs_stored
  attrs$Call <- c(attrs_stored$Call, Call)
  attrs$seed <- c(attrs$seed, seed)
  attributes(aus) <- c(dim=list(dim(aus)), attrs[setdiff(names(attrs),"dim")])
  aus
}
