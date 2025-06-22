#' Make a larger CA from a DHF and a smaller CA
#'
#' Function to construct a larger CA by replacing levels of a DHF with columns of a smaller CA
#'
#' @rdname DHF2CA
#'
#' @aliases DHF2CA
#'
#' @usage DHF2CA(P, D, v = max(D) + 1, unique=TRUE, ...)
#'
#' @param P an M x k DHF of strength t with w levels and at least v classes
#'          (a PHF of strength t with w levels is a special case)
#' @param D an N x w CA of strength t with \code{v} levels
#' @param v number of levels of \code{D},
#'          automatically detected in case of coding as 0,...,\code{v}-
#' @param unique logical; if \code{FALSE}, does not remove non-unique rows from the outcome
#' @param ... further arguments to function \code{\link{maxconstant}}
#'
#' @section Details:
#' A PHF is a perfect hash family, a DHF a distributed hash family (weaker condition).\cr
#' Function \code{DHF2CA} yields a CA(chi + (N-rho)*M, t, k, v) from replacing the levels
#' of the DHF (or the PHF) with the columns of the CA, leaving out the constant rows of the CA
#' and adding a few (chi) back in, if needed, where chi=max(0, rho - (M - 1)(v − rho)),
#' according to Corollary 2.4 of Colbourn and Torres-Jimenez (2010).\cr
#' The function uses function \code{maxconstant} to maximize the constant rows of the CA.
#'
#' The function does not check for the strengths of the ingoing arrays. Hence, it is the users
#' responsibility to ensure that these are as needed.
#'
#' @returns a CA(chi + (N-rho) * M, t, k, v) (matrix, not of class \code{ca}, as there are no
#' guarantees regarding the strength, because guarding this is left to the user or to calling functions).
#'
#' @references Colbourn, Doughterty and Horsley (2019) and references therein, as well as Colbourn and Torres-Jimenez 2010.
#'
#' @examples
#' #########################################
#' ## DHF2CA
#' #########################################
#'
#' ## example from Table 1 of Colbourn and Torres-Jimenez
#' P <- t(SCA_Bose(11)[,1:3])
#' dim(P) ## levels 0, ..., 10
#' D <- paleyCA(3,11) ## 2-level, one constant row
#' dim(D)
#' (chi <- max(0, 1 - 2*(2-1)))
#' aus <- DHF2CA(P, D)
#' dim(aus)  ## 0 + 3*(12-1) rows, 121 columns
#' eCAN(3, 121, 2)  ## optimal
#'
#' # further examples are demos only, they are not for good use cases
#'
#' # a DHF[4,10,4,4,2], i.e., 4 rows, 10 columns, 4 levels, strength 4, 2 partitions
#' # from Colbourn et al. 2021, Table 1
#' P <- cbind(rep(1,4), c(1, rep(2,3)), c(1, rep(3,3)),
#'            rbind(2, rbind(1:3, c(2:3,1), c(3, 1:2))),
#'            rbind(3, rbind(1:3, c(3, 1:2), c(2:3, 1))),
#'            rep(4,4))
#' ## CA with v=2 levels, 4 columns, strength 4 (--> full factorial)
#' D <- as.matrix(expand.grid(0:1,0:1,0:1,0:1))    ## two rows can be made constant
#' (chi <- max(0, 2 - 3*(2-2)))
#' aus <- DHF2CA(P, D, unique=FALSE)
#' dim(aus)     ## 2 + 4*(16-2)= 58 rows, 10 columns
#' aus <- DHF2CA(P, D)
#' dim(aus)     ## 52 x 10
#'              ## 6 duplicate rows removed
#'              ## (this happens for poorly chosen ingredients)
#' coverage(aus, 4)  ## works
#'
#' # a PHF(4,4,19,9) from the Dwyer catalogue
#' P <- rbind(c(4, 2, 3, 7, 4, 0, 1, 2, 5, 3, 4, 5, 7, 6, 1, 0, 1, 8, 8),
#'            c(6, 2, 8, 4, 2, 1, 1, 0, 7, 5, 3, 8, 7, 4, 0, 4, 5, 6, 8),
#'            c(7, 1, 4, 0, 2, 5, 8, 5, 1, 7, 5, 3, 6, 4, 0, 3, 2, 6, 8),
#'            c(4, 8, 5, 0, 0, 1, 2, 5, 1, 1, 7, 7, 2, 4, 3, 8, 6, 3, 0))
#' # we need a strength 4 CA with nine columns
#' D <- paleyCA(4, 9)   ## is current best, has two constant rows
#' aus <- DHF2CA(P, D, unique=FALSE)
#' (chi <- max(0, 2 - (4-1)*(2-2)))
#' dim(aus)     ## not competitive, although D is optimal 2 + (24-2)*4
#' aus <- DHF2CA(P, D)
#' dim(aus)     ## 38 duplicate rows removed, still relevantly worse than current best
#' Ns(4, 19, 2)
#'
#' ## best for 6-level is in DWYERcat, 2768 runs for 9 columns in 6 levels
#' \dontrun{
#'   pathGH <- "https://raw.githubusercontent.com/aadwyer/CA_Database/main/Repository/CA"
#'   # the instruction line starts with a comment character and is ignored,
#'   # therefore ninstruct=0
#'   D <- readCA(paste0(pathGH, "/CA_2768_4_9_6.txt"), ignore.chars=c("]", "["),
#'               ninstruct=0, sep=",")
#'   ## DHF2CA needs quite a while for this case
#'   ## activated one_is_enough argument for maxconstant for run time reasons
#'   aus <- DHF2CA(P, D, one_is_enough=TRUE)  ##
#'   (chi <- max(0, 1 - (2-1)*(6-1)))
#'   dim(aus)  ## 4*2767 = 11068, much more than the 5422 in the Colbourn tables
#'   coverage(aus, 4)
#' }
#'

#' @export
DHF2CA <- function(P, D, v=max(D)+1, unique=TRUE, ...){
  ## function for use by other functions
  ## skips all the checks, except for matrix
  stopifnot(is.matrix(P), is.matrix(D))
  stopifnot(is.numeric(P), is.numeric(D))
  w <- levels.no(P)[1] ## assuming that all rows have the same number of levels
  ##                 but that does not really matter, as long as they are meant to refer
  ##                 to the same CA
  if (min(P)==0) P <- P + 1  ## assuming that levels are either 0 to w-1 or 1 to w
  M <- nrow(P); N <- nrow(D)
  k <- ncol(P)
  levs <- sort(unique(c(D)))
  if (min(levs) > 0){
    ## make levels start with 0
    D <- D - min(levs)
    levs <- levs - min(levs)
  }
  Drem <- maxconstant(D, verbose=2, ...)
  nconst <- length(attr(Drem, "constant_rows")$row_set_list)

  ## nconst is the rho of Corollary 2.4 of Colbourn and Torres-Jimenez
  chi <- max(0, nconst - (M-1)*(v-nconst))
      ## if chi > 0, some constant rows have to be added back in

## Yr the levels that are not constant in the CA used for row r
## YMp1 the levels that were constant for all rows (there are such levels for chi>0)

  const <- c(Drem[1:nconst,1])
  Y <- setdiff(levs, const)
  Drem <- Drem[-(1:nconst),]  ## matrix to use in every step
  aus <- matrix(NA, (N-nconst)*M, k)
  for (i in 1:M){
    ## the constant row(s) are modified from row to row
    aus[(i-1)*(N-nconst) + 1:(N-nconst),] <- (Drem[,P[i,]] + (i - 1)*nconst) %% v
    if (chi>0)
      Y <- union(Y, setdiff(levs, (const + (i-1)*nconst)%%v))
            ## const are the originally constant levels of D
            ## Y holds all levels that were already nonconstant
  }
  if (chi>0){
    YMp1 <- setdiff(levs, Y)
    if (!length(YMp1)==chi) warning("wrong length for YMp1")
    aus <- do.call(rbind, c(as.list(YMp1), list(aus)))
  }
  ## unique is useful for poor choices of P and/or D
  if (unique) return(unique(aus)) else return(aus)
}
