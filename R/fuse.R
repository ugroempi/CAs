#' function to reduce the number of levels by 1
#'
#' reduces the number of levels by 1,
#' achieving a reduction of the number of rows by 2.
#' fuse can also reduce the number of levels by more than 1,
#' reducing the number of rows by 2 for each reduction of 1 level.
#'
#' @rdname fuse
#'
#' @aliases fuse
#'
#' @usage fuse(D, vprior, vpost, ...)
#'
#' @param D a CA with \code{vprior} levels in each column (integers, starting a 0 or 1)
#' @param vprior the number of levels prior to applying fuse
#' @param vpost the number of levels after applying fuse, must be smaller than vprior
#' @param ... currently not used
#'
#' @returns an array of the same coverage properties and the same number of columns
#' as \code{D} and at most \code{nrow(D) - 2 * (vprior - vpost)} rows
#' with \code{vpost} levels per column
#'
#' @section Details:
#' Note that there is a variant of fusion, implemented in function \code{\link{fuseBose}},
#' for which the number of runs can be reduced by three; this is applicable for strength 2
#' with \code{vprior} a prime or prime power and up to \code{vprior + 1} columns in \code{D}.
#'
#' @seealso [fuseBose()] for the case mentioned in the Details section
#'
#' @references Colbourn et al. (2010, CKRS; Lemma 3.1)
#'
#' @examples
#' set.seed(1212) ## for array randomization
#' D <- lhs::createBush(11, 12)
#' table(D)
#' dim(D)
#' D10 <- fuse(D, 11, 10) ## reduce by two rows and one level
#' table(D10)
#' dim(D10)
#' coverage(D, 3)
#' coverage(D10, 3)
#'
swapvals <- function(M, c, val1, val2){
  ## c can be a single integer
  ## or a vector of integers
  hilf <- M[,c]
  stopifnot(all(c(val1,val2) %in% hilf))
  eins <- which(hilf==val1)
  zwei <- which(hilf==val2)
  hilf[eins] <- val2
  hilf[zwei] <- val1
  M[,c] <- hilf
  M
}

#'@export
fuse <- function(D, vprior, vpost, ...){
  stopifnot(is.numeric(vprior))
  stopifnot(is.numeric(vpost))
  if (!vpost < vprior) stop("function fuse reduces number of levels, \n vprior <= vpost is inadequate")
  if (!is.matrix(D)) D <- as.matrix(D)
  stopifnot(is.numeric(D))  ## D must have levels 0 to vprior-1 or 1 to vprior
  k <- ncol(D)
  lev <- levels.no(D)
  if (length(table(lev))>1) stop("function fuse works for symmetric CAs only")
  nlev <- lev[1]
  if (!nlev==vprior) stop("D does not have vprior levels per column")
  for (j in 1:(vprior-vpost)){
    maxlev <- max(D)
    ## how can it be made sure that the 3 run reduction is achieved
    ## for strength 2 with k < maxlev and maxlev-1 a prime power?
    ## according to CKRS Lemma 3.1
    for (c in 1:k){
      if (!D[1,c]==maxlev)
      D <- swapvals(D, c, D[1,c], maxlev) ## first row has maximum number of levels
    }
    D <- D[-1,]
    ## eliminate next row
    elim2 <- D[1,]
    D <- D[-1,]
    repl <- which(elim2==maxlev)  ## arbitrary values for these columns
    for (c in 1:k){
      posmax <- which(D[,c]==maxlev)
      npos <- length(posmax)
      start0 <- maxlev < lev[1]  ## values 0 to maxlev-1 (TRUE) or not
      if (c %in% repl) D[posmax,c] <-
        rep((1-start0):(maxlev-1),npos%/%(nlev-1)+1)[1:npos]
      else D[posmax,c] <- elim2[c]
    }
  }
  D
}
