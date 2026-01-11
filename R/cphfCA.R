## code for creating a CA from the CPHF
#' Making CA from CPHF
#'
#' according to the process outlined in Wagner, Colbourn and Simos (2022)
#'
#' @rdname cphfCA
#'
#' @aliases cphfCA
#' @aliases N_cphfCA
#' @aliases k_cphfCA
#' @aliases WCS
#'
#' @usage cphfCA(t, k, v, start0 = TRUE, ...)
#' @usage N_cphfCA(t, k, v)
#' @usage k_cphfCA(t, N, v)
#' @usage WCS(cphf, q, t, type="2022", ...)
#'
#' @param t integer: target strength
#' @param k integer: requested number of columns
#' @param v integer: number of levels, must be a prime power
#' @param start0 logical: should the levels of the output array start with 0 (otherwise with 1)?
#' @param N integer: affordable run size
#' @param cphf matrix with integer entries that can be expanded into permutation vectors of length \code{t}; must be a CPHF
#' @param q integer: number of levels, must be a prime power
#' @param type character string providing the type of construction; for the built-in CPHFs, \code{"2022"} is the correct entry.
#'        It handles CPHFs for which the first entry of the permutation vector refers to the highest power,
#'        and which were constructed with Galois fields from SageMath or the Python \code{numbthy} module
#'        (relevant only for prime powers 8, 9, 16, 25; also implemented in this package via the internal function \code{mygf})\cr
#'        Any type different from the default requests use of the built-in Galois fields of R package \code{\link{lhs}}.
#' @param ... currently not used
#'
#' @returns Function \code{cphfCA} returns a matrix of class \code{ca} with attributes,
#' which is a covering array of strength 3, 4, 5 or 6.\cr
#' Function \code{WCS} returns a matrix with entries in 0,...,v-1, which is a covering
#' array of strength \code{t}, if a valid \code{cphf} for strength \code{t} is provided;
#' the function is not primarily meant for direct use,
#' but is the workhorse behind \code{cphfCA}.
#'
#' Functions \code{N_cphfCA} and \code{k_cphfCA} return the required run size or the achievable number of columns, respectively.
#'
#' @section Details:
#' Function \code{cphfCA} implements CAs from the Wagner et al. (2022) construction
#' for length \code{t} h-tuples to create length \code{v^t} permutation
#' vectors for \code{v} levels and strength \code{t}.
#' The workhorse function \code{WCS} uses a stored covering perfect hash family (CPHF); at present,
#' only the CPHFs from Wagner, Colbourn and Simos (2022) are implemented via the function
#' \code{cphfCA} (see \code{\link{WCS_CPHFs}}
#' and \code{\link{CPHFcat}}).\cr
#' Users who have their own valid CPHF (not SCPHF)
#' that is not part of this package can use the function for creating a CA. Note that
#' the creation of a CA can take a long time for large cases; for example, creation of a
#' CA(501121, 4, 665, 17) from the relevant CPHF with 6 rows (\code{6*(17^4-1)+1=501121})
#' took approximately 140 minutes on the authors Windows computer with 32GB RAM
#' (which appeared to be the limiting factor). This package is not optimized to handle cases
#' as large as this.
#'
#' The construction relies on a Galois field GF(\code{v}).
#' The CPHFs from Wagner et al. (2022) for \code{v=8}, \code{v=9}, \code{v=16} and \code{v=25}
#' require a different Galois field than otherwise used in this package; these deviating
#' Galois fields are created with the internal function \code{mygf}; function \code{WCS} takes
#' the argument \code{type}, with which these Galois fields are activated (per default), but can
#' also be deactivated in favor of the Galois fields of package \code{\link{lhs}}
#' by specifying any type different from \code{"2022"}.\cr
#' The data frame \code{\link{CPHFcat}} holds the overview information on which constructions
#' are implemented for function \code{cphfCA}.
#'
#' @references Wagner, Colbourn and Simos (2022)
#'
#' @examples
#' # applying the construction for small CPHF from the Wagner et al. (2022) paper
#'
#' # strength 3
#' head(CPHFcat[which(CPHFcat$v==4),])
#' ## the CPHF for the first row
#' WCS_CPHFs[["3"]][["4"]][["71"]]
#' D <- cphfCA(3,71,4)
#' dim(D)
#' \dontrun{
#'   coverage(D,3)
#' }
#' Ns(3,71,4)   ## size is current best but not eCAN
#'
#' ## the CPHF for the 91st row has flexible entries
#' CPHFcat[91,]
#' rowSums(is.na(WCS_CPHFs[["4"]][["3"]][["3504"]]))
#' ## many flexible values in the last row,
#' ## eCAN claims one row less for this construction
#' \dontrun{
#'   ## creation of this CA is slow because of the many columns
#'   D <- cphfCA(4,3504,3)  ## about 3 minutes on author's machine
#'   dim(D)
#' }
#' Ns(4,3504,3)   ## by far the best implemented
#'
#' # strength 6
#' Ns(6, 10, 3)  ## by far best implemented,
#'               ## slightly better array from group who does not share their arrays
#' dim(D <- bestCA(6,10,3))
#' coverage(D, 6)
#'
#'
#' @importFrom sfsmisc digitsBase

#' @export
cphfCA <- function(t, k, v, start0=TRUE, ...){
  Call <- sys.call()
  stopifnot(is.numeric(t), is.numeric(k), is.numeric(v))
  zeilen <- CPHFcat[CPHFcat$t>=t & CPHFcat$k>=k & CPHFcat$v==v,]
  ## lower strength is at the top,
  ## higher strength at the bottom
  zrows <- nrow(zeilen)
  if (zrows == 0) stop("The request cannot be served by function cphfCA.")
  zeile <- zeilen[zrows + 1 - which.min(rev(zeilen$N)), ,drop=FALSE]
  tc <- as.character(zeile$t)
  kc <- as.character(zeile$k)
  vc <- as.character(zeile$v)
  P <- WCS_CPHFs[[tc]][[vc]][[kc]]
  ## arbitrary fix NA values at v
  if (zeile$hasNA) P[which(is.na(P))] <- v
    aus <- WCS(P, zeile$v, zeile$t)
  if (!start0) aus <- aus + 1
  aus <- aus[,1:k]
  class(aus) <- c("ca", class(aus))
  attr(aus, "Call") <- Call
  attr(aus, "origin") <- "Walker Colbourn Simos 2022 CPHF"
  attr(aus, "t") <- zeile$t
  aus
}

#' @export
WCS <- function(cphf, q, t, type="2022", ...){
  # cphf a strength t CPHF, with tuples of length t
  # in q^t levels
  ## matrix with the beta vectors in the order t-1 top
  ## to 0 bottom
  ## type="2022" is the Wagner Colbourn Simos paper of 2022

  ## powers 8, 9, 16, 25 are from Python/SAGEmath (like in SMC paper for 8 and 9)
  if (q %in% c(8,9,16,25) && type=="2022") gf <- mygf(q, type="2006") else
    gf <- lhs::create_galois_field(q)

  ## coefficients for low powers at the top
  ## nothing to be appended, as no SCPHF

  ## lspace provides the h tuples for the CA columns
  ## separately per cphf row
  ## top is for the highest power, bottom for the lowest
  ## the list of txk matrices of permutation vectors
  lspace <- lapply(1:nrow(cphf), function(obj)
          sfsmisc::digitsBase(cphf[obj,], q, ndigits=t))
  ## betaspace contains the beta coefficients
  ## top/left is for the highest power, bottom for the lowest,
  ##      which is the default of digitsBase
  betaspace <- t(sfsmisc::digitsBase(0:(q^t-1), q, ndigits=t))
  unique(do.call("rbind",
                 lapply(lspace, function(obj)
                   gf_matmult(betaspace, obj, gf))))
}

#' @export
N_cphfCA <- function(t, k, v){
  stopifnot(is.numeric(t), is.numeric(k), is.numeric(v))
  zeilen <- CPHFcat[CPHFcat$t==t & CPHFcat$k>=k & CPHFcat$v==v,]
  ## lower strength is at the top,
  ## higher strength at the bottom
  zrows <- nrow(zeilen)
  if (zrows == 0) return(NA)
  zeile <- zeilen[zrows + 1 - which.min(rev(zeilen$N)), ,drop=FALSE]
  min(zeilen$N)
}

#' @export
k_cphfCA  <- function(t, N, v){
  hilf <- CPHFcat[CPHFcat$t==t & CPHFcat$N<=N &
                               CPHFcat$v==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  max(hilf[,"k"])
}

