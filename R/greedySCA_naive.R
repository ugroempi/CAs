#' Greedy random creation of an SCA
#'
#' for k columns and strength t
#'
#' @rdname greedySCA_naive
#'
#' @aliases greedySCA_naive
#' @aliases iter_greedySCA_naive
#'
#' @usage greedySCA_naive(k, t, seed=NULL, postopt=TRUE, ...)
#' @usage iter_greedySCA_naive(niter, k, t, seed_init=NULL, verbose=TRUE, ...)
#'
#' @param k positive integer (at least t); the number of columns
#' @param t the strength for which the coverage of all permutations of length t must be guaranteed for all \code{combn(k,t)} t-element subsets of 1,...,k
#' @param seed an optional positive integer seed for random sampling; specifying a seed makes the array generation reproducible
#' @param postopt logical; default TRUE, i.e., rows with all elements redundant are identified and removed. If this takes too long, change to FALSE
#' @param niter number of independent attempts
#' @param seed_init initial seed for iterative application
#' @param verbose logical (default TRUE), controlling how much printed output is created
#' @param ... currently not used
#'
#' @returns a strength t SCA with k columns; all rows contain the integers 1 to k.
#'
#' @section Details:
#' Function \code{greedySCA_naive} starts from t! random permutations of the the
#' integers 1,...,k and adds a further random permutation,
#' until all orders of all t-tuples are covered. For \code{postopt=TRUE},
#' it then identifies the redundant elements and removes all rows that
#' contain redundant elements only.
#'
#' Function \code{iter_greedySCA_naive} makes \code{niter} independent attempts and
#' improves the best of those by removing rows that contain redundant elements only.
#'
#' @section Use of AI:
#' Claude 4 was involved in the development of these functions.
#'
#' @seealso [greedySCA_Kuhn()] for the greedy method by Kuhn et al. (2012) and [greedySCA_TJ()] for the greedy stage method of Torres-Jimenez et al. (2022)
#'
#' @author Ulrike Groemping
#'
#' @examples
#' ## without post-optimization
#' nrow(greedySCA_naive(5, 3, seed=2323, postopt=FALSE))
#' ## relatively unlucky
#' nrow(greedySCA_naive(5, 3, seed=2415, postopt=FALSE))
#' ## relatively lucky
#' nrow(A <- greedySCA_naive(5, 3, seed=9456, postopt=FALSE))
#' coverageSCA(A, 3)
#' coverageSCA(A, 4)
#'
#' ## with post-optimization
#' ## differences become less dramatic (but still relevant)
#' nrow(greedySCA_naive(5, 3, seed=2323))
#' ## relatively unlucky
#' nrow(greedySCA_naive(5, 3, seed=2415))
#' ## relatively lucky
#' nrow(A <- greedySCA_naive(5, 3, seed=9456))
#' coverageSCA(A, 3)
#' coverageSCA(A, 4) ## worse coverage of 4-sequences than for larger SCA
#'
#' ## iterative with unlucky start
#' Aiter <- iter_greedySCA_naive(10, 5, 3, seed_init=2415)
#' nrow(Aiter)
#' ## niter=20 better than niter=10
#' Aiter <- iter_greedySCA_naive(20, 5, 3, seed_init=2415)
#' nrow(Aiter)
#'
#' ## further improvement is possible by more substantial attempt
#' nrow(reduce_rows_iterative_complete(Aiter, 3))

#' @export
greedySCA_naive <- function(k, t, seed=NULL, postopt=TRUE, ...){
  stopifnot(is.numeric(k), is.numeric(t))
  k <- round(k); t <- round(t)
  stopifnot(length(k)==1, length(t)==1)
  stopifnot(k>=t, t > 1)
  stopifnot(is.logical(postopt))
  combs <- combn(k, t)
  perms <- combinat::permn(t)
  if (!is.null(seed)) set.seed(seed)

  ## All t-tuples to cover
  all_subsequences <- do.call(rbind, unlist(lapply(1:ncol(combs), function(obj)
    lapply(perms, function(obj2) combs[,obj][obj2])), recursive = FALSE))
  checkbits <- rep(FALSE, nrow(all_subsequences))
  cur <- t(sapply(1:factorial(t), function(obj) sample(k)))

  for (i in 1:nrow(all_subsequences)){
    seq <- all_subsequences[i,]
    for (j in 1:nrow(cur)){
      if (funcheck(cur[j,], seq)){
        checkbits[i] <- TRUE
        break
      }
    }
  }

  while(!all(checkbits)){
    cur <- rbind(cur, sample(k))
    #print(dim(cur))
    #print(cur[nrow(cur),])
    for (i in which(!checkbits)){
      seq <- all_subsequences[i,]
      if (funcheck(cur[nrow(cur),], seq))
        checkbits[i] <- TRUE
    }
  }
  if (postopt){
    hilf <- identify_needed(cur, t)
    cur <- cur[!rowSums(hilf)==0,]
  }
  cur
}

#' @export
iter_greedySCA_naive <- function(niter, k, t, seed_init=NULL, verbose=TRUE, ...){
  stopifnot(is.numeric(niter), is.numeric(k), is.numeric(t))
  niter <- round(niter); k <- round(k); t <- round(t)
  stopifnot(length(niter)==1, length(k)==1, length(t)==1)
  stopifnot(niter>0, k>=t, t>1)
  nzeilen <- Inf
  curmat <- NULL
  for (i in 1:niter){
    if (verbose) cat("=== Iteration: ", i, "\n")
    if (i==1) seed <- seed_init else seed <- NULL
    hilf <- greedySCA_naive(k, t, postopt=FALSE, seed=seed)
    if (nrow(hilf) < nzeilen){
      curmat <- hilf
      nzeilen <- nrow(curmat)
      if (verbose) cat("=== Current no. of rows:", nzeilen, "rows\n")
    }
  }
  hilf <- identify_needed(curmat, t)
  if (verbose) cat("===", length(which(rowSums(hilf)==0)), "redundant rows removed\n")
  curmat[!rowSums(hilf)==0,]
}
