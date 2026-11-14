#' function for optimal 2-level CAs of strength t with t+2 columns and for optimal
#' v-level CAs of strength t with t+1 columns in v^t runs
#'
#' implements the construction by Johnson and Entringer (19xx) and
#' the Zero Sum construction reported by Colbourn (20xx)
#'
#' @rdname JE
#' @aliases JE
#' @aliases k_JE
#' @aliases N_JE
#' @aliases ZeroSum
#' @aliases k_ZeroSum
#' @aliases N_ZeroSum
#'
#' @usage JE(t=NULL, k=NULL)
#' @usage k_JE(t=NULL, N=NULL)
#' @usage N_JE(t=NULL, k=NULL)
#' @usage ZeroSum(t, k, v)
#' @usage k_ZeroSum(t, N, v)
#' @usage N_ZeroSum(t, k, v)
#'
#' @param t requested strength
#' @param k number of columns
#' @param N number of runs
#'
#' @details
#' Functions \code{ZeroSum} and \code{JE} return the
#' optimum CAs in strength \code{t} with
#' \code{k=t+1} columns (zero sum construction, valid for any \code{v})
#' or \code{k=t+2} columns for \code{v=2} (Johnson-Entringer construction).
#' As they are optimal for exactly these settings and better arrays
#' can be found for smaller values of \code{k}, the \code{N_}-functions
#' return \code{NA} for \code{k} less than the respective maximum, in order to provide useful input for the \code{\link{Ns}}
#' and \code{\link{ks}} functions.
#'
#' The design creation functions increase
#' \code{t} if \code{k} is larger than possible for the stated \code{t}, with a message.
#' They return \code{k} columns only if more were possible. In such cases,
#' the construction is (of course) not optimal and also neve the
#' best possible implemented construction.
#'
#' @return \code{JE} returns a matrix of class \code{ca} with
#' \code{k} columns (at most \code{t+2}) in 2 levels each
#' and \code{2^t} runs in each column, with coverage strength \code{t}.
#'
#' \code{N_JE} and \code{k_JE} return the \code{N} for a given \code{t}
#' and \code{k} or the \code{k} for a given \code{t} and \code{N},
#' respectively.
#'
#' \code{ZeroSum} returns a matrix of class \code{ca} with
#' \code{k} columns (at most \code{t+1}) in \code{v} levels each
#' and \code{v^t} runs in each column, with coverage strength \code{t}.
#'
#' \code{N_ZeroSum} and \code{k_ZeroSum} return the \code{N} for a given
#' \code{t}, \code{v} and \code{k} or the \code{k} for a given \code{t},
#' \code{v} and \code{N} respectively.
#'
#' @references Colbourn and Dinitz (2007, III 7.16 for Zero Sum),
#' Johnson and Entringer (1989)
#'
#' @examples
#' D <- JE(k=8)
#' N_JE(k=8)
#' dim(D)
#' k_JE(t=3,N=66)
#' N_JE(t=3,k=6)
#' N_JE(t=3, k=5)
#' dim(JE(t=3,k=5))

#' @export
N_JE <- function(t=NULL, k=NULL) {
  if (is.null(t) && is.null(k))
    stop("At least one of k and t must be specified")
  if (is.null(t)){
    t <- k-2
    return(floor(2^(t+2)/3))
  } else ## t was not NULL
    { if (is.null(k) || k == t+2) return(floor(2^(t+2)/3)) else
      return(NA)
      }
}

#' @export
k_JE <- function(t=NULL, N=NULL) {
  if (is.null(t) && is.null(N)) stop("At least one of N and t must be specified")
  if (is.null(t))
    t <- floor(log2(ceiling(N*3/4))) else{ ## t was not NULL
      if (!is.null(N) && N < floor(2^(t+2)/3)) return(NA)
    }
    ## now t is defined, and no invalid combination
    return(t+2)
}

#' @export
JE <- function(t=NULL, k=NULL){
  Call <- sys.call()
  ## k is the number of columns
  ## t is the strength

  # Input validation with meaningful error messages
  if (is.null(k) && is.null(t))
    stop("At least one of k and t must be specified")
  ## if one is NULL, it is calculated from the other
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k != as.integer(k))
    stop("k must be a single positive integer")

  if (!is.null(t)){
    if (!is.numeric(t) || length(t) != 1 || t < 2 || t != as.integer(t))
      stop("t must be a single integer >= 2")
  }
  if (!is.null(k) && !is.null(t)){
    if (k > t+2) {
      t <- k-2 ## make stronger for reaching k
      message(paste0("The strength was increased, \nbecause there is no JE array with the requested strength and ", k, " columns."))
    }
  }
  if (is.null(k)) k <- t+2
  if (is.null(t)) t <- k-2
  hilf <- ffCA(t+2, 2)
  ham <- rowSums(hilf)
  posmin <- which.min(table(ham %% 3))
  aus <- hilf[which(ham%%3==posmin-1),1:k]
  class(aus) <- c("ca", class(aus))
  attr(aus, "t") <- t
  attr(aus, "Call") <- Call
  attr(aus, "origin") <- "Johnson and Entringer"
  attr(aus, "comment") <- "optimal for k=t+2"
  aus
}

#' @export
N_ZeroSum <- function(t, k, v) {
  if (!is.numeric(t) || length(t) != 1 || t < 1 || t != as.integer(t))
    stop("t must be a single positive integer")
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k != as.integer(k))
    stop("k must be a single positive integer")
  if (!is.numeric(v) || length(v) != 1 || v < 1 || v != as.integer(v))
    stop("v must be a single positive integer")
  if (!k == t+1) return(NA)
  return(v^t)
}

#' @export
k_ZeroSum <- function(t, N, v) {
  if (!is.numeric(t) || length(t) != 1 || t < 1 || t != as.integer(t))
    stop("t must be a single positive integer")
  if (!is.numeric(N) || length(N) != 1 || N < 1 || N != as.integer(N))
    stop("N must be a single positive integer")
  if (!is.numeric(v) || length(v) != 1 || v < 1 || v != as.integer(v))
    stop("v must be a single positive integer")
  if (N < v^t) return(NA)
  if (N >= v^(t+1)) message(paste0("a larger t would be possible with N=", N, ", \nthe returned value refers to the stated ", t))
  t+1
}

#' @export
ZeroSum <- function(t, k, v){
  Call <- sys.call()
  ## k is the number of columns
  ## t is the strength

  # Input validation with meaningful error messages
  if (!is.numeric(t) || length(t) != 1 || t < 1 || t != as.integer(t))
    stop("t must be a single positive integer")
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k != as.integer(k))
    stop("k must be a single positive integer")
  if (!is.numeric(v) || length(v) != 1 || v < 1 || v != as.integer(v))
    stop("v must be a single positive integer")

  if (!is.null(t) && !is.null(k))
    if (k > t+1) {
      t <- k-1
      message(paste0("The strength was increased, \nbecause there is no ZeroSum array with the requested strength and ", k, " columns."))
    }
  hilf <- ffCA(t,v)
  aus <- cbind(hilf, (-rowSums(hilf))%%v)
  class(aus) <- c("ca", class(aus))
  attr(aus, "t") <- t
  attr(aus, "Call") <- Call
  attr(aus, "origin") <- "ZeroSum"
  attr(aus, "comment") <- "optimal for k=t+1"
  aus
}

