#' Making mixed-level CAs
#'
#' either from a given CA D, by makes unused levels of some columns flexible (
#' arbitrary strength), or for strength 2 only,
#' from required numbers of levels, by a construction generalized from Sherwood (2008).
#' In either case, the result is treated by function postopNCK for reducing its
#' run size.
#'
#' @rdname MCA
#'
#' @aliases MCA2
#' @aliases CA_to_MCA
#'
#' @usage MCA2(nlevels, D=NULL, outerRetry=10, innerRetry=1, seed=NULL, ...)
#' @usage CA_to_MCA(D, cs, tolevs, t=attr(D, "t"), outerRetry=10, innerRetry=1, seed=NULL, ...)
#'
#' @param nlevels a vector of numbers of levels, or a data.frame that is the table
#'       of that vector, with levels in decreasing order and column names \code{level}
#'       and \code{frequency}
#' @param D for \code{MCA2}, a CA to be used for expanding; if \code{NULL}, it is
#' automatically determined; it is meant for experts only, and it may have mixed levels\cr
#' for \code{CA_to_MCA}, a strength \code{t} CA (it may also have mixed levels)
#' @param outerRetry passed to \code{\link{postopNCK}}; if \code{0},
#'       run size reduction will not be attempted.
#' @param innerRetry positive integer value, passed to \code{\link{postopNCK}};
#'       it is sometimes but not always beneficial to increase it versus the default 1
#' @param seed the seed for making \code{\link{postopNCK}} reproducible;\cr
#'        if \code{NULL}, it will be randomly determined and reported with the result
#' @param D a strength \code{t} CA (it may have mixed levels)
#' @param cs integer vector of column numbers
#' @param tolevs integer vector of target numbers of levels for the columns indicated
#'        by \code{cs}; must not be larger than the number of levels of the related
#'        columns in \code{cs}, and must not be smaller than 2
#' @param t integer, the strength of \code{D} (smaller strength would be permitted,
#'        but would not make sense)
#' @param ... currently not used
#'
#' @section Details:
#' Function \code{MCA2} implements Groemping's (2025) generalization of the first
#' construction of Sherwood (2008) for creating mixed covering arrays (MCAs).
#' The algorithm is described in Groemping (2025). Its key idea is to
#' start from a uniform CA, make superfluous levels flexible, and expand
#' columns with enough flexible entries into more columns of the same number of levels,
#' using a recursive substitution of flexible values with ordered designs.
#' The result is then post-processed by function \code{\link{postopNCK}},
#' unless this is suppressed by the argument \code{outerRetry=0}. Post-processing
#' can also be done separately by directly using function \code{\link{postopNCK}}
#' on the result of function \code{MCA2}, regardless whether it was created with
#' or without post-processing.
#'
#' Function \code{CA_to_MCA} sets unused levels of the columns in \code{cs} to \code{NA}, i.e.,
#' makes them flexible, as, e.g., described in Moura et al. (2003). Subsequently, the function
#' \code{\link{postopNCK}} applies the method of Nayery et al. (2013) for
#' making entire rows flexible and removing them. Repeated runs can remove more rows.
#' The \code{outerRetry} argument gives the number of calls to \code{\link{postopNCK}}.
#' if the bottom rows of the outcome still look promising, the function can be called
#' again for removing even more rows.
#'
#' Without specifying a seed, the results will differ from call to call.
#' The seed, even
#' if not explicitly specified, is stored with the output, so that the result will
#' always be reproducible, if desired.
#'
#' The run size optimization can take a long time. It can be suppressed by setting
#' \code{outerRetry} to zero, and this should be done for large settings.
#' Per default, the optimization is switched on. Its progress is reported
#' by interim messages. The functions react beneficial to user interrupts: They
#' return the result from the previous successful \emph{outer} retry (make sure the outer
#' retry has finished before escaping calculations!).
#'
#' @returns Both functions return a matrix of class \code{ca} with attributes of the ingoing matrix (dimensions modified),
#'    and attributes \code{Call} and \code{seed} modified, and with flexible values in most columns.
#'
#' @section Warning:
#' There is not much experience yet with the post-optimization performance.
#' The defaults for the related parameters may change in the future.
#'
#' @references Groemping (2025), Moura et al. (2003), Nayeri et al. (2013), Sherwood (2008)
#'
#' @author Ulrike Groemping
#'
#' @examples
#' ##################################################
#' ### MCA2
#' ##################################################
#' ## small example without run size optimization
#' ## 14 runs (one better than CAgen for this setting)
#' D <- MCA2(c(rep(2,30), rep(3,2)), outerRetry=0)
#' attributes(D)
#' coverage(D, 2)
#'
#' \dontrun{
#' ## optimize run size (using a good seed)
#' D <- MCA2(c(rep(2,30), rep(3,2)), seed=17369, outerRetry=3)
#' coverage(D, 2)
#' attributes(D)
#' }
#' ## reasonably quickly reduces to 11 runs on my Windows machine
#' ## lower bound is 9 runs (3*3)
#'
#' ## create a large design without trying a run size optimization
#' D <- MCA2(c(rep(2,70), rep(3,80), rep(4,8), rep(6,2)),
#'          outerRetry=0)
#' dim(D)
#' coverage(D, 2)
#' ## alternative call for the same setting
#' D <- MCA2(data.frame(levels=c(6,4,3,2), frequency=c(2,16,80,70)),
#'          outerRetry=0)
#' \dontrun{
#' ## this runs for a long while, and it is not known whether a reduction
#' ## can be achieved; the code can be interrupted and returns the latest
#' ## improved array (if any)
#' Doptimized <- postopNCK(D, 2)
#' }
#'
#' ## calls with more extreme numbers of levels
#' D <- MCA2(c(13,11,rep(4,8), rep(2,20)), outerRetry=0)
#' ## 169 x 30
#' ## minimum run size is known to be 13*11=143
#' ## obtained relatively fast with postopNCK
#'
#' D <- MCA2(c(13,rep(4,8), rep(3,12), rep(2,20)), outerRetry=0)
#' ## 169 x 41
#' ## lower bound for run size is 13 * 4 = 52
#' ## postopNCK with seed 3318, innerRetry=3:
#' ##    169 quickly shrunk to 64 runs (after 9 outer retries)
#' ##    reasonably quickly to 59 runs (after 13 outer retries),
#' ##    then slowly to 53 runs (after 19 outer retries)
#' ##    stopping short of reaching the known lower bound 52
#' ##    within the default settings for outer retries
#'
#' ## use the argument D
#' D <- MCA2(c(rep(5,1), rep(4,2), rep(3,7), rep(2,20)),
#'          outerRetry=0)
#'
#' ########################################################
#' ### CA_to_MCA
#' ########################################################
#' # a small example
#' ## six 3-level columns in 33 runs
#' D <- bestCA(3,6,3)
#' ## to three 3-level and three 2-level columns in 27 runs
#' Dmixed <- CA_to_MCA(D, cs=4:6, tolev=rep(2, 3), t=3, outerRetry=0)
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
CA_to_MCA <- function(D, cs, tolevs, t=attr(D, "t"),
                      outerRetry=10, innerRetry=1, seed=NULL, ...){
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
    message("D already has the optimum number of rows.")
    attrs <- attrs_stored
    attrs$Call <- c(attrs_stored$Call, Call)
    attributes(D) <- attrs
    return(D)
  }
  Dnow <- D
  if (outerRetry >=1)
  tryCatch({
      aus <- postopNCK(Dnow, t, outerRetry = outerRetry,
                       innerRetry = innerRetry, seed=seed)
      latest_result <- aus
      Dnow <- aus
  },
  interrupt = function(e) {
    attr(latest_result, "rowOrder") <- NULL
    attrs <- attrs_stored
    attrs$Call <- c(attrs_stored$Call, Call)
    attrs$seed <- c(attrs$seed, seed)
    attributes(latest_result) <- c(dim=list(dim(latest_result)),
                                   attrs[setdiff(names(attrs),"dim")])
    return(latest_result)
  })
  else aus <- D
  if (nrow(aus) == N){
    if (outerRetry > 0) message("The number of rows of D could not be reduced.")
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

#' @export
MCA2 <- function(nlevels, D=NULL, outerRetry=10, innerRetry=1,
                 seed=NULL, ...){
  Call <- sys.call()
  if (is.null(seed)) seed <- sample(32000, 1)
  if (!is.list(nlevels)){
    if (!is.numeric(nlevels))
      stop("if not a list or a data.frame, \nnlevels must be a vector of numbers of levels")
  if (is.numeric(nlevels)){
    k <- length(nlevels)
    tab <- table(nlevels)
    if (!(k==sum(tab))) stop("nlevels must have non-NA entries only")
    nlevels <- data.frame(
      levels=as.numeric(rev(names(tab))),
      frequency=rev(unname(as.numeric(tab))))
    }
  } ## now nlevels is a list or data.frame
  ## correct format
  stopifnot(c("levels", "frequency") %in% names(nlevels))
  ## decreasing order
  if (!all(nlevels$levels==sort(nlevels$levels, decreasing=TRUE)))
    stop("levels in nlevels must be sorted in decreasing order")
  ## unique row elements
  if (!nrow(nlevels)==length(unique(nlevels$levels)))
    stop("nlevels$levels must not have duplicate entries")
  ## optimum N
  optN <- ifelse(nlevels$frequency[1]==2, nlevels$levels[1^2],
                 prod(nlevels$levels[1:2]))

  ## find a suitable construction
  hilf <- find_mix_construction_separate(nlevels)
  ## create the ingoing array
  if (is.null(D))
      D <- bestCA(2, hilf$k, hilf$v)
  attrs <- attributes(D)
  ## apply the construction
  ## outerseparate re-sorts the levels of D in the order of increasing frequency
  ## in order to allow as many columns as possible for the expansion
  aus <- outerseparate(D, nlevels, ...)
  latest_result <- aus
  message("found array with ", nrow(aus), " rows")

  ## try to improve the run size by postopNCK
  ## made early-stop friendly
  Dnow <- aus
  if (outerRetry > 0){
  message("working on improving the run size")
  message("on interrupt, the current best result is returned")
  message("there will be at most ", outerRetry, " outer loop steps")
  tryCatch({
      aus <- postopNCK(Dnow, 2, innerRetry = innerRetry, outerRetry=outerRetry, seed=seed)
      latest_result <- aus
      Dnow <- aus
    },
    interrupt=function(e){
      attrnow <- attributes(latest_result)
      namesvorher <- names(attrs)
      namesnow <- names(attrnow)
      attrs$dim <- attrnow$dim
      attrs$dimnames <- attrnow$dimnames
      attrs[setdiff(names(attrnow), names(attrs))] <- attrnow[setdiff(names(attrnow), names(attrs))]
      attributes(latest_result) <- attrs
      attr(latest_result, "Call") <- c(Call, attr(latest_result, "Call"))
      attr(latest_result, "eCAN") <- NULL
      attr(latest_result, "PCAstatus") <- NULL
      attr(latest_result, "rowOrder") <- NULL  ## no longer relevant
      return(latest_result)
  })
  }
  attrnow <- attributes(aus)
  namesvorher <- names(attrs)
  namesnow <- names(attrnow)
  attrs$dim <- attrnow$dim
  attrs$dimnames <- attrnow$dimnames
  attrs[setdiff(names(attrnow), names(attrs))] <- attrnow[setdiff(names(attrnow), names(attrs))]
  attributes(aus) <- attrs
  attr(aus, "Call") <- c(Call, attr(aus, "Call"))
  attr(aus, "eCAN") <- NULL  ## no longer relevant
  attr(aus, "rowOrder") <- NULL  ## no longer relevant
  attr(aus, "PCAstatus") <- NULL  ## will be destroyed
  return(aus)
}
