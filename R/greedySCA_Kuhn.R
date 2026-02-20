#' Greedy random creation of an SCA
#'
#' for k columns and strength t according to an algorithm by Kuhn et al. (2012)
#'
#' @rdname greedySCA_Kuhn
#'
#' @aliases greedySCA_Kuhn
#' @aliases iter_greedySCA_Kuhn
#'
#' @usage greedySCA_Kuhn(k, t, nsamp=10, postopt=TRUE, seed=NULL, ...)
#' @usage iter_greedySCA_Kuhn(niter, k, t, nsamp=10, postopt=FALSE,
#'     seed_init=NULL, verbose=TRUE, ...)
#'
#' @param k positive integer (at least t); the number of columns
#' @param t the strength for which the coverage of all permutations of length t must be guaranteed for all \code{combn(k,t)} t-element subsets of 1,...,k
#' @param nsamp the number of runs sampled in each step of the algorithm (the N of the paper); the default 10 will likely be too small for realistic applications
#' @param postopt logical; default depends on function; if TRUE, the function \code{\link{reduce_rows_iterative_complete}} tries to make rows redundant and remove them. If this takes too long, \code{postopt} should be FALSE.
#' @param seed an optional positive integer seed for random sampling; specifying a seed makes the array generation reproducible
#' @param niter the number of independent attempts to make
#' @param seed_init an optional positive integer seed for the first iteration; specifying a seed makes the array generation reproducible
#' @param verbose logical (default TRUE), controlling how much printed output is created
#' @param ... currently not used
#'
#' @returns a strength t SCA with k columns; all rows contain the integers 1 to k.
#'
#' @section Details:
#' The function \code{greedySCA_Kuhn} implements the construction of Kuhn et al. (2012). It starts
#'   from t! random permutations of the the integers 1,...,k and adds a further random permutation,
#'   until all orders of all t-tuples are covered.\cr
#'   For \code{postopt=TRUE}, it then *identifies the redundant elements and removes all rows that
#'   contain redundant elements only (would be easiest, but at present not done)* applies iterative run size reduction
#'   using function \code{reduce_rows_iterative_complete}.
#'
#' The function \code{iter_greedySCA_Kuhn} iteratively repeats the process and keeps the best
#' outcome.
#'
#' @section Use of AI:
#' Claude 4 was involved in the development of these functions (mainly for debugging).
#'
#' @seealso [greedySCA_naive()] for a naive random assembly of an SCA and [greedySCA_TJ()] for the graph-based greedy stage method of Torres-Jimenez et al. (2022)
#'
#' @author Ulrike Groemping
#'
#' @references Kuhn et al. (2012)
#'
#' @examples
#' ## without post-optimization
#' nrow(greedySCA_Kuhn(5, 3, seed=2323, postopt=FALSE))
#' ## relatively unlucky
#' nrow(greedySCA_Kuhn(5, 3, seed=2415, postopt=FALSE))
#' ## relatively lucky
#' nrow(A <- greedySCA_Kuhn(5, 3, seed=16, postopt=FALSE))
#' coverageSCA(A, 3)
#' coverageSCA(A, 4)
#'
#' ## with post-optimization
#' ## differences become smaller (still relevant esp. for larger cases)
#' nrow(greedySCA_Kuhn(5, 3, seed=2323))
#' ## relatively unlucky
#' nrow(greedySCA_Kuhn(5, 3, seed=2415))
#' ## relatively lucky
#' nrow(A <- greedySCA_Kuhn(5, 3, seed=16))
#' coverageSCA(A, 3)
#' coverageSCA(A, 4) ## worse coverage of 4-sequences than for larger SCA
#'
#' ## iteration
#' nrow(iter_greedySCA_Kuhn(10, 5, 3, seed_init=2415)) ## start unlucky
#'   ## larger niter brings the run size down to 8

#' @export
greedySCA_Kuhn <- function(k, t, nsamp=10, postopt=TRUE, seed=NULL, ...){
  ## work with numbers 1 to k (not 0 to k-1)
  call <- sys.call()
  combs <- combn(k,t)
  # perms <- arrangements::permutations(n=t)
  # perms <- lapply(1:factorial(t), function(obj) perms[obj,])
  ## ??? Is the order of permutations relevant for this method?
  perms <- combinat::permn(t)  ## a list
  aus <- do.call(rbind, unlist(lapply(1:ncol(combs), function(obj)
    lapply(perms, function(obj2) combs[,obj][obj2])), recursive=FALSE))
  checkbits <- rep(0, nrow(aus))
  ## initialize with first two rows (consistently more successful for small test cases)
  ts <- rbind(1:k, k:1)  ## initial order
  initiallycovereds <- apply(ts, 1,
                             function(obj) which(apply(aus, 1, function(obj2) funcheck(obj,obj2))))
  checkbits[initiallycovereds] <- 1
  # ts <- matrix(NA, 0, k)
  nuncovered <- sum(1-checkbits)
  ## set a seed for sampling
  ## and for picking a single run from a non-unique best set
  if (is.null(seed)) seed <- sample(32000, 1)
  set.seed(seed)
  zaehl <- 0
  while(nuncovered > 0){
      tc <- replicate(nsamp, sample(k))  ## matrix with nsamp columns
    # else
    #   tc <- t(arrangements::permutations(n=k))  ## matrix with nsamp columns

    nuncovered <- sum(1-checkbits)
    whichuncovereds <- which(checkbits==0)
    newlycovereds <- apply(tc, 2,
                           ### pick the rows from aus[whichuncovereds,] for which
                           ### this column of tc covers the sequence
                           ### obj is a column of tc
                           function(obj) whichuncovereds[which(apply(aus[whichuncovereds,,drop=FALSE],
                                                                     ## obj2 is a row of aus[whichuncovereds,]
                                                                     1, function(obj2) funcheck(obj, obj2)))])
    ## pick the most successful ones
    laengen <- lengths(newlycovereds)
    laengen[is.na(laengen)] <- 0

    if (all(laengen <= 0)) next

    pick <- which.max(laengen)
    #sample(which(laengen==max(laengen)), 1)
    test1 <- tc[,pick]
    ts <- rbind(ts, test1)
    checkbits[newlycovereds[[pick]]] <- 1
    nuncovered <- sum(1-checkbits)
    if (any(checkbits==0)){
      ## add reverse order
      # test2 <- rev(test1);
      ## add inverted order (coincides with reverse for first pair)
      test2 <- k - test1 + 1
      whichuncovereds <- which(checkbits==0)
      newlycovereds <- whichuncovereds[which(apply(aus[whichuncovereds,,drop=FALSE],
                                                   1, function(obj2) funcheck(test2,obj2)))]
      ts <- rbind(ts, test2)
      checkbits[newlycovereds] <- 1
      nuncovered <- sum(1-checkbits)
    }}
  if (postopt){
    print(paste("Number of rows before postopt: ", nrow(ts)))
    ## this is more demanding than just removing redundant rows
    ## might downgrade to do only that ?
    ts <- reduce_rows_iterative_complete(ts, t)
    print(paste("Number of rows after postopt: ", nrow(ts)))
  }
  rownames(ts) <- NULL
  class(ts) <- c("sca", class(ts))
  attr(ts, "call") <- call
  attr(ts, "seed") <- seed
  ts
}

#' @export
iter_greedySCA_Kuhn <- function(niter, k, t, nsamp=10, postopt=FALSE,
                                seed_init=NULL, verbose=TRUE, ...){
  nzeilen <- Inf
  curmat <- NULL
  for (i in 1:niter){
    if (verbose) cat("=== Iteration: ", i, "\n")
    if (i==1) seed <- seed_init else seed <- NULL
    hilf <- greedySCA_Kuhn(k, t, nsamp=nsamp, postopt=postopt, seed=seed)
    if (nrow(hilf) < nzeilen){
      curmat <- hilf
      nzeilen <- nrow(curmat)
      if (verbose) cat("=== Current no. of rows:", nzeilen, "rows\n")
    }
  }
  curmat
}
