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
#' @aliases N_upper_MCA
#' @aliases projBoseMCA
#'
#' @usage MCA2(nlevels, D=NULL, outerRetry=10, innerRetry=1, seed=NULL, ...)
#' @usage CA_to_MCA(D, cs, tolevs, t=attr(D, "t"), outerRetry=10, innerRetry=1, seed=NULL, ...)
#' @usage N_upper_MCA(nlevels, t=2, internet = TRUE, ...)
#' @usage projBoseMCA(nlevels, t=2, ...)
#'
#' @param nlevels a vector of numbers of levels, or a data.frame that is the table
#'       of that vector, with levels in decreasing order and column names \code{level}
#'       and \code{frequency}; function \code{MCA2} always sorts the numbers of levels
#'       into decreeasing order
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
#'        but would not make sense) or the requested strength
#' @param internet logical: is an internet connection available (for DWYER or NIST download)
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
#' Function \code{N_upper_MCA} reports an upper bound on the run number for an MCA
#' of the setting. That bound is obtained by using function \code{\link{bestN}} for
#' the maximum number of levels with the total number of columns. Depending on the situation,
#' the bound may be sharp (e.g., for eight five-level columns with two 3-level columns)
#' or dramatically too large (e.g., for one eight-level column with ten 2-level columns).
#'
#' Function \code{projBoseMCA} implements Theorem 2.3 of Colbourn (2008) as well as
#' Corollary 2.2 of the same paper, which also references Stevens, Ling and Mendelsohn
#' (2002, their Theorem 2.6). It works for two different numbers of levels only.
#' It is most powerful for q+1 columns at q-1 levels with a single column at q+1 levels,
#' where q is a prime power (the Corollary 2.2 situation);
#' in these cases the run size is optimal ((q-1)*(q+1)=q^2-1).
#'
#' @returns Functions \code{MCA2} and \code{CA_to_MCA} return a matrix of class \code{ca}
#'    with attributes of the ingoing matrix (dimensions modified),
#'    and attributes \code{Call} and \code{seed} modified,
#'    and often with flexible values in most columns. For function \code{MCA},
#'    the number of levels is in decreasing order
#'    (regardless of the order specified in \code{nlevels}).\cr
#'    Function \code{N_upper_MCA} returns an upper bound for the run size of an MCA
#'    obtainable via the constructions (named integer, obtained by \code{\link{bestN}},
#'    where the size can be obtained by using the CA of the construction indicated
#'    by the name as the \code{D} in function \code{CA_to_MCA}).\cr
#'    Function \code{projBoseMCA} returns a matrix of class \code{ca} with attributes;
#'    for numeric \code{nlevels}, the columns are in the order of \code{nlevels}, otherwise
#'    the columns with more levels come before the columns with fewer levels.
#'
#' @section Warning:
#' There is not much experience yet with the post-optimization performance.
#' The defaults for the related parameters may change in the future.
#'
#' @references Colbourn (2008), Groemping (2025), Moura et al. (2003), Nayeri et al. (2013),
#' Sherwood (2008), Stevens, Ling and Mendelsohn (2002)
#'
#' @author Ulrike Groemping
#'
#' @examples
#' ##################################################
#' ### MCA2
#' ##################################################
#' ## small example with fast optimization to global optimum
#' ## which is 2*8
#' ## 16 runs
#' D <- MCA2(c(rep(2,10), 8))
#' dim(D)
#' coverage(D, 2)
#' ## N_upper_MCA is very pessimistic
#' N_upper_MCA(c(rep(2,10), 8))
#'
#' ## another case for which N_upper_MCA seems to be tight
#' N_upper_MCA(c(rep(5,8), 3,3))
#' D <- CA_to_MCA(CS_LCDST(10,5), 9:10, c(3,3), outerRetry=0)
#' dim(D)
#' head(D); tail(D)
#'
#' \dontrun{
#' ## optimize run size (using a good seed)
#' D <- MCA2(c(rep(2,30), rep(3,2)), seed=17369, outerRetry=3)
#' coverage(D, 2)
#' attributes(D)
#' ## reasonably quickly reduces to 11 runs on my Windows machine
#' ## lower bound is 9 runs (3*3)
#' ## CAgen yields 14 runs with IPOG-F
#' }
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
#' ## use the argument D;
#' ## even this naive version, though starting large, is quickly reduced
#' ## to adequate size (sometimes even faster than smaller one)
#' D <- MCA2(c(5, rep(4,2), rep(3,7), rep(2,20)), D=bestCA(2,30,5),
#'          outerRetry=0)
#' \dontrun{
#'    Doptimized <- postopNCK(D, 2)
#' }
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
#' ######################################################
#' ####           projBoseMCA               #############
#' ######################################################
#'
#' ## corollary 2.2 cases (a single number of levels is different,
#' ## and two larger than the others)
#' D <- projBoseMCA(c(rep(4,6),6))
#' head(D)
#' dim(D) ## 24 runs like IPOG-F for CAgen
#'
#' D <- projBoseMCA(c(rep(6,8),8))
#' head(D)
#' dim(D) ## 48 runs, better than the 60 for IPOG-F
#'        ## which is best for CAgen and the 58 of JMP Pro
#'        ## as of 11 August 2025
#'
#' ## the single different number can also be smaller
#' ## then not guaranteed to be optimal
#' D <- projBoseMCA(c(rep(6,8),7))
#' head(D)
#' dim(D) ## 48 runs, better than the 60 for IPOG-F
#'        ##  42 would be a guaranteed optimal
#'
#' D <- projBoseMCA(c(rep(6,8),4))
#' head(D)
#' dim(D) ## 48 runs, better than the 55 for IPOG-F
#'        ##  36 would be a guaranteed optimal
#'
#' ## cases with two groups of alphabet sizes
#' ## the larger alphabet must be more frequent
#' D <- projBoseMCA(c(rep(4,5), rep(3,4)))
#' dim(D)  ## based on q=8 and c=4,
#'         ## there could be up to 9 4-level columns
#' D <- projBoseMCA(c(rep(4,9), rep(3,4)))
#' dim(D)  ## 60 runs, larger than 30 for IPOG-F2
#'
#' ## better use CA_to_MCA for this setting
#' D <- CA_to_MCA(bestCA(2,13,4), 10:13, rep(3,4), outerRetry=0)
#' dim(D)  ## better than CAgen, even without postoptimization
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
  optN <- ifelse(nlevels$frequency[1]==2, nlevels$levels[1]^2,
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

#' @export
N_upper_MCA <- function(nlevels, t=2, internet=TRUE, ...){
  if (!is.list(nlevels))
    if (!is.numeric(nlevels))
      stop("if not a list or a data.frame, \nnlevels must be a vector of numbers of levels")

    if (is.numeric(nlevels)){
      k <- length(nlevels)
      v <- max(nlevels)
    }else{
      k <- sum(nlevels$frequency)
      v <- max(nlevels$levels)
    }
  bestN(t, k, v, internet=internet, ...)
}

#' @export
projBoseMCA <- function(nlevels, t=2, ...){
  if (!is.list(nlevels))
    if (!is.numeric(nlevels))
      stop("if not a list or a data.frame, \nnlevels must be a vector of numbers of levels")

  if (is.numeric(nlevels)){
    nlevels_ordered <- nlevels
    levels <- sort(unique(nlevels), decreasing = TRUE)
    if (length(levels) < 2) stop("uniform CA requested, use projBoseCA")
    if (length(levels) > 2) stop("projBoseMCA not applicable for more than two different numbers of levels")
    nlevels <- data.frame(levels=levels, frequency=sapply(levels, function(obj) sum(nlevels==obj)))
  }else{
    nlevels_ordered <- rep(nlevels$levels, nlevels$frequency)
    levels <- nlevels$levels
    if (length(levels) > 2) stop("projBoseMCA not applicable for more than two different numbers of levels")
    if (!all(sort(levels, decreasing = TRUE)==levels)) stop("data.frame nlevels must be sorted with levels in descending order")
  }
  ## now, nlevels is a data.frame with levels sorted in decreasing order
  c_temp <- min(nlevels$frequency)
     ## c at least c_temp, perhaps more
  s <- nlevels$levels[which(nlevels$frequency==c_temp)]
  levother <- setdiff(nlevels$levels, s)
  if (s > levother && c_temp > 1) stop("combination of c and s not permissible")
  if (s > levother){
    ## the special case of Colbourn 2008 Corollary 2.2
    c <- 1
    if (!(s - levother <= 2)) stop("this setting is not covered by projBoseMCA")
    if (!(levother + 1 %in% primedat$q)) stop("this setting is not covered by projBoseMCA")
    q <- levother + 1
    nqminus1 <- max(nlevels$frequency)
    if (nqminus1 > q + 1) stop("two many columns with ", q, "-1 levels for projBoseMCA")
    aus <- projectionBose(q, c, s)[,c(1:nqminus1, q+2)]
    ## arrange columns in requested order for originally numeric
    ## nlevels
    pos_s <- which.max(nlevels_ordered)
    if (!pos_s==length(nlevels_ordered)){
      if (pos_s==1) aus <- aus[,c(ncol(aus),1:(ncol(aus)-1))]
      else aus <- aus[,c(1:(pos_s-1), ncol(aus), pos_s:(ncol(aus)-1))]
    }
    ## assign attributes
    class(aus) <- c("ca", class(aus))
    attr(aus, "origin") <- "Colbourn 2008 Corollary 2.2 (and Stevens, Ling, Mendelssohn 2002)"
    attr(aus, "t") <- 2
    return(aus)
  }
  ## now s < levother
  ## c_temp = frequency of s
  c <- NA
  if (levother + c_temp %in% primedat$q &&
      max(nlevels$frequency) <= levother + c_temp + 1){
      c <- c_temp
      q <- levother + c
    }else{
      ## now c_temp must be increased for getting
      ## both a prime and enough columns
      repeat{
        c_temp <- c_temp + 1
        if (levother + c_temp %in% primedat$q &&
            max(nlevels$frequency) <= levother + c_temp + 1){
          c <- c_temp
          q <- levother + c
          break
        }
      }
    }
  ## as s < levother
  nqminusc <- max(nlevels$frequency)
  n_s <- min(nlevels$frequency)
  ## create the array with maximum number of columns
  hilf <- projectionBose(q, c, s)[,c(1:nqminusc, (q+2):(q+1+n_s))]
  ## arrange columns in requested order for originally numeric
  ## nlevels
  pos_s <- which(nlevels_ordered==min(nlevels$levels))
  pos_qminusc <- which(nlevels_ordered==max(nlevels$levels))
  aus <- matrix(NA, nrow=nrow(hilf), ncol=ncol(hilf))
  aus[,pos_s] <- hilf[,(nqminusc+1):(ncol(hilf))]
  aus[,pos_qminusc] <- hilf[,1:nqminusc]

  ## assign attributes
  class(aus) <- c("ca", class(aus))
  attr(aus, "origin") <- "Colbourn 2008 thm 2.3"
  attr(aus, "t") <- 2
  aus
}
