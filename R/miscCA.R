#' Incorporate miscellaneous arrays
#'
#' incorporates miscellaneous arrays that are created in the package
#' or obtained from package DoE.base
#'
#' @rdname miscCA
#'
#' @aliases miscCA
#'
#' @usage miscCA(t, k, v, fixNA=TRUE, seed=NULL, maxconstant=FALSE, makePCA=FALSE, ...)
#'
#' @param t integer: the strength
#' @param k integer: the number of columns
#' @param v integer: the number of levels for each column
#' @param fixNA logical: should flexible values be fixed?\cr
#'       If the array has flexible values, \code{fixNA=TRUE}
#'       randomly assigns a fixed value to each flexible value.\cr
#'       Thus, a seed is needed for reproducibility.
#' @param seed \code{NULL} or an integer seed for making the
#'       fixing of NA values reproducible; if \code{seed=NULL},
#'       a seed is randomly obtained and stored in the attribute
#'       \code{fixNA_seed} of the returned object.
#' @param maxconstant logical: should constant rows be maximized ?\cr
#'       ignored for internal arrays, applied for arrays loaded
#'       from package \pkg{\link{DoE.base}} only;
#'       per default, the array is not modified and may not even
#'       have a single constant row
#' @param makePCA logical: should a PCA structure be enforced ?\cr
#'       ignored for internal arrays, applied for arrays loaded
#'       from package \pkg{\link{DoE.base}} only\cr
#'       If both \code{makePCA} and \code{maxconstant} are TRUE,
#'       \code{maxconstant} takes precedence.
#' @param ... currently not used
#'
#' @returns an array of class \code{ca} with attributes
#'
#' @section Details:
#' The function relies on the data frame \code{\link{miscCAcat}} for a
#' description of the available CAs.
#'
#' @references Ji and Yin (2010), Colbourn and Torres-Jimenez (2013), Kuhfeld (w/o year),
#' \pkg{\link{DoE.base}}
#'
#' @examples
#' # a CA from DoE.base
#' miscCA(3, 5, 4)
#'
#' ## without making rows constant and without enforcing PCA
#' attributes(miscCA(3, 5, 4, maxconstant=FALSE))
#'
#' # a CA that is internally stored
#' # it has a single flexible value, which is per default fixed
#' attributes(D <- miscCA(2, 14, 4))
#' Ns(2,14,4)  ## best possible
#'

#' @export
miscCA <- function(t, k, v, fixNA=TRUE, seed=NULL, maxconstant=FALSE, makePCA=FALSE, ...){
  Call <- sys.call()
  hilf <- miscCAcat[miscCAcat[,"t"]==t &
                      miscCAcat[,"k"]>=k &
                      miscCAcat[,"v"]==v,,drop=FALSE]
  if (nrow(hilf)==0) stop("no suitable array is available")
  hilf <- hilf[which.min(hilf$N),,drop=FALSE]
  if (length(grep("DoE.base", hilf$fns))==1){
      aus <- eval(parse(text=hilf$fns))
      origin <- attr(aus, "origin")
      if (makePCA && maxconstant && hilf$nconst < hilf$v){
        message("maxconstant takes precedence over makePCA, both are not possible")
        message("Set maxconstant to FALSE in order to obtain a PCA structure with k1 = ",
        hilf$PCAstatus)
      }
      if (maxconstant){
        if (hilf$nconst==1) aus <- maxconstant(aus, one_is_enough = TRUE) else
        aus <- maxconstant(aus)
        }
      class(aus) <- c("ca", class(aus))
      dimnames(aus) <- NULL
      attr(aus, "t") <- hilf$t
      attr(aus, "origin") <- origin
      if (makePCA && !maxconstant){
        aus <- CA_to_PCA(aus)
        aus <- CA_to_PCA(aus, tryhard = TRUE)
      }
      comment <- ""
      if (!maxconstant) comment <- paste0(hilf$nconst, " constant rows are possible")
      if ((maxconstant || !makePCA) && !hilf$nconst == hilf$v)
        comment <- paste0("PCA structure with k1=", hilf$PCAstatus, " can be achieved")
      if (!comment=="") attr(aus, "comment") <- c(attr(aus, "comment"), comment)
  }else{
    if (length(grep("CohenSA", hilf$fns))==1) aus <- eval(parse(text=hilf$fns))
        else aus <- get(hilf$fns)
      dimnames(aus) <- NULL
      ## it is assumed that the attributes are already suitably set
      ## maxconstant and makePCA are ignored
    }
    ## from DoE.base and internal treated together
    if (ncol(aus) > k){
        attrs <- attributes(aus)
        attrs$dim[2] <- k
        aus <- aus[,1:k]
        attributes(aus) <- attrs
      }
    attr(aus, "Call") <- Call
    ## handle missing (=flexible) values
    if (hilf$hasNA && any(is.na(aus))){
        ## the "any" part is there because it might be the case
        ## that the column(s) with NA were removed because of k < ncol
      if (!fixNA)
        attr(aus, "flexible") <- list(value=NA, profile=colSums(is.na(aus)))
      else{
        if (is.null(seed)) seed <- sample(1:32000, 1)
        attr(aus, "fixNA_seed") <- seed
        tobefilled <- which(is.na(aus))
        nfill <- length(tobefilled)
        set.seed(seed)
        aus[tobefilled] <- sample(0:(v-1), nfill, replace=(nfill > v))
      }
    }
  aus
}
