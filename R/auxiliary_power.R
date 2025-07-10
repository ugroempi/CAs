#' Auxiliary functions for Colbourn Torres-Jimenez Power construction
#'
#' These functions support the implementation of the power constructions
#' by Colbourn and Torres-Jimenez
#'
#' @rdname auxpower
#'
#' @aliases createDHF
#' @aliases Tred
#'
#' @usage createDHF(oa, toa, tdhf, vdhf=tdhf)
#' @usage Tred(DHF, rs, reduce)
#'
#' @param oa   an orthogonal array
#' @param toa  strength of the orthogonal array
#' @param tdhf requested strength of the DHF
#' @param vdhf number of classes of the DHF, the default yields a PHF
#' @param DHF a DHF or DHHF, for which the number of levels for
#'            rows indexed by elements of \code{rs}
#'            is to be reduced by \code{reduce}
#' @param rs   vector of row indices whose number of levels is to be reduced
#' @param reduce number of levels to reduce
#'
#' @returns Function \code{createDHF} creates a DHF or PHF of minimum
#' number of rows for the combination of \code{toa}, \code{tdhf} and
#' \code{vdhf}.\cr
#' Function \code{Tred} reduces the number of levels for the specified rows
#' by the amount stated in \code{reduce}, aiming for retaining the
#' maximum number of columns.
#'
#' @section Details:
#' Function \code{createDHF} creates a distributed hash family according to
#' Lemma 1.1 or Lemma 1.2 of Colbourn and Torres-Jimenez (2010).
#'
#' Function \code{Tred} applies Lemma 3.1 of Colbourn and Torres-Jimenez (2010).
#'
#' @references Colbourn and Torres-Jimenez (2010)
#'
#' @examples
#' ## a DHF of strength 3 with 2 classes
#' ## from the Bush array of strength 2
#' ## (which coincides with the Bose array)
#' mydhf <- createDHF(SCA_Busht(11, 2), 2, 3, 2)
#' dim(mydhf)
#' mydhf[,1:6]
#' mydhf[,116:121]
#' ## reduce the third row by 6 levels to 5 levels
#' mydhhf <- Tred(mydhf, 3, 6)
#' dim(mydhhf)
#' mydhhf[,50:55]
#'
#' ## a PHF of strength 4 from a strength 2 OA
#' myphf19 <- createDHF(SCA_Bose(19), 2, 4)
#' dim(myphf19)
#'
#' ## a DHF of strength 4 with only 2 classes
#' ## from a strength 2 OA
#' mydhf19 <- createDHF(SCA_Bose(19), 2, 4, 2)
#' dim(mydhf19)
#'

#' @export
createDHF <- function(oa, toa, tdhf, vdhf=tdhf){
  ## if vdhf is not specified, a phf is generated
  k <- ncol(oa)
  bound <- (toa - 1)*Turan(tdhf, vdhf)
  stopifnot(k > bound)
  return(t(oa[,1:(bound+1)]))
}

#' @export
Tred <- function(DHF, rs, reduce){
  start0 <- TRUE
  if (min(DHF)==1){
    start0 <- FALSE
    DHF <- DHF-1
  }
  vs <- levels.no(t(DHF))
  if (!length(unique(vs[rs]))==1) warning("rs contains rows with different levels.")
  for (r in rs){
    freqs <- table(DHF[r,])
    freqs <- sort(freqs)
    levs <- as.numeric(names(freqs))
    vprior <- max(levs) + 1
    ## remove least frequent levels, and rearrange levels to be from 0 to max
    DHF <- DHF[, !DHF[r,] %in% levs[1:reduce]]
    ## remove gaps
    for (i in 0:(vprior - reduce - 1)){
       if (!i %in% DHF[r,]){
         subtract <- min(DHF[r, DHF[r,] > i]) - i
         DHF[r, DHF[r,] > i] <- DHF[r, DHF[r,] > i] - subtract
       }
    }
  }
  DHF + as.numeric(!start0)
}
