## code for creating a CA from the CPHF
#' Making CA from SCPHF
#'
#' according to the process outlined in Sherwood, Martirosyan and Colbourn (2006)
#' (non-reduced permutation vectors).
#'
#' @rdname scphfCA
#'
#' @aliases scphfCA
#' @aliases SMC
#' @aliases N_scphfCA
#' @aliases k_scphfCA
#'
#' @usage scphfCA(t, k, v, start0 = TRUE, ...)
#' @usage SMC(cphf, q, t, type="2006", ...)
#' @usage N_scphfCA(t, k, v)
#' @usage k_scphfCA(t, N, v)
#'
#' @param t integer: target strength
#' @param k integer: requested number of columns
#' @param v integer: number of levels, must be a prime power
#' @param start0 logical: should the levels of the output array start with 0 (otherwise with 1)?
#' @param cphf matrix, which is a covering perfect hash family (CPHF)
#' @param q integer: number of levels, prime power
#' @param type string: "2006" refers to the original implementation by Sherwood,
#'     Martirosyan and Colbourn (2006), which uses non-default Galois fields for q=8,9;\cr
#'     string "2009" refers to SCPHFs from Walker and Colbourn (2009) which also uses non-default Galois fields for q=8,9, and uses a different one for \code{type="2006"} for q=9;\cr
#'     anything else uses the built-in Galois fields of package \code{\link{lhs}}
#' @param N integer: affordable run size
#' @param ... currently not used
#'
#' @returns Function \code{scphfCA} returns a matrix of class \code{ca} with attributes,
#' which is a covering array of strength 3 or 4.\cr
#' The function \code{SMC} returns a matrix with entries in 0,...,q-1, which is a covering
#' array of strength \code{t}, if a valid \code{cphf} for strength \code{t} is provided;
#' the function is not primarily meant for direct use,
#' but is the workhorse behind \code{scphfCA}. It can be used with user-provided CPHFs.
#'
#' Functions \code{N_scphfCA} and \code{k_scphfCA} return the required run size or the achievable number of columns, respectively.
#'
#' @section Details:
#' Function \code{scphfCA} implements CAs from the Sherwood et al. (2006) construction
#' for length \code{t-1} h-tuples to create length \code{v^t} permutation
#' vectors for \code{v} levels and strength \code{t}.
#' The workhorse function \code{SMC} uses a stored covering perfect hash family of a
#' specific form (SCPHF). The built-in SCPHFs are from Sherwood et al. (2006)
#' and from Lanus as provided in the Dwyer (2024) GitHub repository (see \code{\link{SMC_SCPHFs}}
#' and \code{\link{CL_SCPHFs}}). Users who want to use their own SCPHF can implement it with the function
#' \code{SMC}; large cases can take a long time to construct.\cr
#' The construction relies on a Galois field GF(\code{v}).
#' The SCPHFs from Sherwood et al. (2006) for \code{v=8} or \code{v=9} require a different Galois field
#' than otherwise used in this package (created with internal function \code{mygf}); the \code{type} argument
#' of function \code{SMC} can be used for requesting those (type="2006") or the default ones (type="2018").\cr
#' The data frame \code{\link{SCPHFcat}} holds the overview information on which constructions
#' are implemented.
#'
#' @references Sherwood, Martirosyan and Colbourn (2006), Walker II and Colbourn (2009), Colbourn and Lanus (2018), Colbourn, Lanus and Sarkar (2018), Dwyer (2024)
#'
#' @examples
#' # applying the construction for small SCPHF from the Sherwood et al. (2006) paper
#'
#' # strength 3
#' head(SCPHFcat[which(SCPHFcat$v==4),])
#' ## the SCPHF for the first row (type 2006)
#' SMC_SCPHFs[["3"]][["4"]][["16"]]
#' D <- scphfCA(3,16,4)
#' dim(D)
#' coverage(D,3)
#' eCAN(3,16,4)   ## size matches eCAN
#'
#' ## the SCPHF for the third row (type 2018)
#' CL_SCPHFs[["3"]][["4"]][["55"]]
#' D <- scphfCA(3,55,4)
#' dim(D)
#' \dontrun{coverage(D,3, parallel=4)}
#' Ns(3,55,4)   ## best implemented
#'
#' # strength 4
#' Ns(4, 10, 3)  ## best possible
#'
#'
#' @importFrom sfsmisc digitsBase

#' @export
scphfCA <- function(t, k, v, start0=TRUE, ...){
  Call <- sys.call()
  stopifnot(is.numeric(t), is.numeric(k), is.numeric(v))
  zeilen <- SCPHFcat[SCPHFcat$t>=t & SCPHFcat$k>=k & SCPHFcat$v==v,]
  ## lower strength is at the top,
  ## higher strength at the bottom
  zrows <- nrow(zeilen)
  if (zrows == 0) stop("The request cannot be served by function scphfCA.")
  zeile <- zeilen[zrows + 1 - which.min(rev(zeilen$N)), ,drop=FALSE]
  vc <- as.character(zeile$v); tc <- as.character(zeile$t); kc <- as.character(zeile$k)
  if (zeile$type=="2006")
      aus <- SMC(SMC_SCPHFs[[tc]][[vc]][[kc]], zeile$v, zeile$t)
  if (zeile$type=="2009")
    aus <- SMC(WC_SCPHFs[[tc]][[vc]][[kc]], zeile$v, zeile$t, type="2009")
  if (zeile$type=="2018")
    aus <- SMC(CL_SCPHFs[[tc]][[vc]][[kc]], zeile$v, zeile$t, type="2018")
  if (!start0) aus <- aus + 1
  aus <- aus[,1:k]
  class(aus) <- c("ca", class(aus))
  attr(aus, "Call") <- Call
  attr(aus, "origin") <- ifelse(zeile$type=="2006",
                                "Sherwood, Martirosyan, Colbourn (2006) SCPHF",
                                ifelse(zeile$type=="2009",
                                       "Walker and Colbourn (2009) SCPHF",
                                       "Lanus and co-authors (2018) related SCPHF"))
  attr(aus, "t") <- t
  aus
}

#' @export
SMC <- function(cphf, q, t, type="2006", ...){
  # cphf a strength t SCPHF, with tuples of length t
  # in q^t levels
  ## matrix with the beta vectors in the order 0 top
  ## to t-1 bottom
  ## type="2006" is the original Sherwood et al. paper
  ## type="2009" is the Walker and Colbourn 2009 paper
  ## type="2018" are the Lanus arrays as downloaded from Dwyer
  ## the type is relevant for the gf for 8 and 9 levels

  if (!q %in% c(8,9,16,25) || !type %in% c("2006", "2009"))
    gf <- lhs::create_galois_field(q)
  else gf <- mygf(q, type=type)

  ## coefficients for low powers at the top
  ## 1 will be prepended

  ## lspace provides the h tuples for the CA columns
  ## separately per cphf row
  lspace <- lapply(1:nrow(cphf), function(obj)
    rbind(1,
          sfsmisc::digitsBase(cphf[obj,], q, ndigits=t-1)[(t-1):1,]))
  ## betaspace contains the beta coefficients
  betaspace <- t(sfsmisc::digitsBase(0:(q^t-1), q, ndigits=t)[t:1,])
  unique(do.call("rbind",
   lapply(lspace, function(obj)
     gf_matmult(betaspace, obj, gf))))
}

#' @export
N_scphfCA <- function(t, k, v){
  stopifnot(is.numeric(t), is.numeric(k), is.numeric(v))
  zeilen <- SCPHFcat[SCPHFcat$t==t & SCPHFcat$k>=k & SCPHFcat$v==v,]
  ## lower strength is at the top,
  ## higher strength at the bottom
  zrows <- nrow(zeilen)
  if (zrows == 0) return(NA)
  zeile <- zeilen[zrows + 1 - which.min(rev(zeilen$N)), ,drop=FALSE]
  min(zeilen$N)
}

#' @export
k_scphfCA  <- function(t, N, v){
  hilf <- SCPHFcat[SCPHFcat$t==t & SCPHFcat$N<=N &
                               SCPHFcat$v==v,,drop=FALSE]
  if (nrow(hilf)==0) return(NA)
  max(hilf[,"k"])
}

