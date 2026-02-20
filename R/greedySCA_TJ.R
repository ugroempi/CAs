#' Greedy random creation of an SCA
#'
#' for k columns and strength t according to an algorithm by Torres-Jimenez et al. (2022)
#'
#' @rdname greedySCA_TJ
#'
#' @aliases greedySCA_TJ
#' @aliases iter_greedySCA_TJ
#'
#' @usage greedySCA_TJ(k, t, postopt=TRUE, seed=NULL, ...)
#' @usage iter_greedySCA_TJ(niter, k, t, postopt=FALSE, seed_init=NULL, verbose=TRUE, ...)
#'
#' @param k positive integer (at least t); the number of columns
#' @param t the strength for which the coverage of all permutations of length t must be guaranteed for all \code{combn(k,t)} t-element subsets of 1,...,k
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
#' The function \code{greedySCA_TJ} implements the construction of Torres-Jimenez et al. (2022). It
#'   uses an approach based on directed acyclic graphs and relies on the R package igraph for that
#'   purpose.\cr
#'   For \code{postopt=TRUE}, it then *identifies the redundant elements and removes all rows that
#'   contain redundant elements only (would be easiest, but at present not done)* applies iterative run size reduction
#'   using function \code{reduce_rows_iterative_complete}.
#'
#' The function \code{iter_greedySCA_TJ} iteratively repeats the process and keeps the best
#' outcome.
#'
#' Torres-Jimenez et al. (2022) describe a 3-stage approach, the third stage of which is
#' simulated annealing. Simulated annealing is invoked if the iterative reduction process is unable
#' to find a row that can be removed without deteriorating coverage. While this stage is needed for
#' achieving the very good performance that the paper reports on, it is often excessively slow at least
#' in the implementation of this package.
#'
#' @section Use of AI:
#' Claude 4 was involved in the development of these functions (mainly for debugging).
#'
#' @seealso [greedySCA_naive()] for a naive greedy random assembly of an SCA and [greedySCA_Kuhn()] for the greedy method by Kuhn et al. (2012)
#'
#' @author Ulrike Groemping
#'
#' @references Torres-Jimenez et al. (2022)
#'
#' @examples
#' ## without post-optimization
#' nrow(greedySCA_TJ(5, 3, seed=2323, postopt=FALSE))
#' ## relatively unlucky
#' nrow(greedySCA_TJ(5, 3, seed=222, postopt=FALSE))
#' ## relatively lucky
#' nrow(A <- greedySCA_TJ(5, 3, seed=8881, postopt=FALSE))
#' coverageSCA(A, 3)
#' coverageSCA(A, 4)
#'
#' ## with post-optimization
#' ## differences become smaller (still relevant esp. for larger cases)
#' nrow(greedySCA_TJ(5, 3, seed=2323))
#' ## relatively unlucky
#' nrow(greedySCA_TJ(5, 3, seed=222))
#' ## very lucky
#' nrow(A <- greedySCA_TJ(5, 3, seed=8881))
#' coverageSCA(A, 3)
#' coverageSCA(A, 4) ## worse coverage of 4-sequences than for larger SCA
#'
#' ## iteration
#' nrow(iter_greedySCA_TJ(10, 5, 3, seed_init=222)) ## start unlucky
#'   ## larger niter brings the run size down to 8

#' @export
greedySCA_TJ <- function(k, t, postopt=TRUE, seed=NULL, ...){
  ## work with numbers 1 to k (not 0 to k-1)
  call <- sys.call()
  combs <- combn(k,t)
  perms <- combinat::permn(t)
  # perms <- arrangements::permutations(n=t)
  # perms <- lapply(1:factorial(t), function(obj) perms[obj,])
  ## the to-be-covered t-tuples
  aus <- do.call(rbind, unlist(lapply(1:ncol(combs), function(obj)
    lapply(perms, function(obj2) combs[,obj][obj2])), recursive=FALSE))
  ## checkbits indicate that not yet covered (1 = uncovered, 0 = covered)
  checkbits <- rep(1, nrow(aus))
  sca <- rbind(1:k, k:1)  ## initial orders
  initiallycovereds <- apply(sca, 1,
                             function(obj) which(apply(aus, 1, function(obj2) funcheck(obj,obj2))))
  checkbits[initiallycovereds] <- 0

  whichuncovereds <- which(checkbits==1)
  hilf <- integer(0)

  if (!is.null(seed)) set.seed(seed)

  while(length(whichuncovereds) > 0){
    ## add missing elements to hilf in random order
    hilf <- c(hilf, sample(setdiff(whichuncovereds, hilf)))

    ## populate a new DAG
    G <- igraph::make_graph(edges=rbind(aus[hilf[1],-t],
                                aus[hilf[1],-1]), n=k)

    if (length(hilf) >= 2) {
      ## extend the graph, as long as no loop is introduced
      for (j in 2:length(hilf)){
        G2 <- igraph::add.edges(G, rbind(aus[hilf[j],-t],
                                 aus[hilf[j],-1]))
        if (!igraph::is_dag(G2))
          break
        else {
          G <- G2
        }
      }
    }

    cand <- igraph::topo_sort(G)

    ## find newly covered tuples
    newlycovereds <- whichuncovereds[which(apply(aus[whichuncovereds,,drop=FALSE],
                                                 1, function(obj2) funcheck(cand,obj2)))]

    if (length(newlycovereds)==0) stop("no improvement")
    ## update checkbits
    checkbits[newlycovereds] <- 0

    ## update whichuncovereds and hilf
    whichuncovereds <- setdiff(whichuncovereds, newlycovereds)
    hilf <- setdiff(hilf, newlycovereds)

    ## add candidate to sca
    sca <- rbind(sca, cand)
  }## end of while loop

  if (postopt) sca <- reduce_rows_iterative_complete(sca, t)
  rownames(sca) <- NULL
  class(sca) <- c("sca", class(sca))
  attr(sca, "call") <- call
  attr(sca, "seed") <- seed
  sca
}

#' @export
iter_greedySCA_TJ <- function(niter, k, t, seed_init=NULL, verbose=TRUE, ...){
  nzeilen <- Inf
  curmat <- NULL
  for (i in 1:niter){
    if (verbose) cat("=== Iteration: ", i, "\n")
    if (i==1) seed <- seed_init else seed <- NULL
    hilf <- greedySCA_TJ(k, t, seed=seed)
    if (nrow(hilf) < nzeilen){
      curmat <- hilf
      nzeilen <- nrow(curmat)
      if (verbose) cat("=== Current no. of rows:", nzeilen, "rows\n")
    }
  }
  curmat
}
