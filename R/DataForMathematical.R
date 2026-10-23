globalVariables(c("CYCLOTOMYcat", "PALEYcat", "rmsCOPIEScat"))

#' Data for function mathematical constructions
#'
#' The data.frames CYCLOTOMYcat and PALEYcat support the identification of valid
#' scenarii for use of the respective constructions, and the actual implementation.
#'
#' @docType data
#'
#' @format
#' \code{CYCLOTOMYcat} is a data.frame in 213 rows
#' with columns \code{t} (strength),
#' \code{k} (number of CA columns), \code{v} (number of levels), \code{N} (number of CA rows),
#' \code{q} (prime or prime power behind the construcion), \code{type} (the construction type),
#' and the code for the workhorse function that the interface function \code{\link{cyclotomyCA}}
#' will use.
#'
#' \code{PALEYcat} is a data.frame in 178 rows with columns \code{t} (strength),
#' \code{k} (number of CA columns), \code{v} (number of levels), \code{N} (number of CA rows),
#' \code{q} (prime or prime power behind the construcion), and \code{nconst} (number of constant
#' rows that can be achieved by function \code{\link{maxconstant}}).
#' Here, \code{v} is always 2. The code is for the workhorse function
#' \code{paley_mat}, if necessary with rows added or removed, will be used by
#' the interface function \code{\link{paleyCA}}.
#'
#' @source
#' The reference for the cyclotomy construction is Colbourn (2010).
#' Some parameters in the paper were found to be wrong. Thus,
#' the parameters in \code{CYCLOTOMYcat} are from the Colbourn tables,
#' where hope is that they are correct.\cr
#' One instance from the Colbourn tables was removed,
#' as its apparent underlying q is not a prime,
#' but 13*229 so that either a different construction was used
#' or the entry carries wrong information. For three further instances,
#' the run size was corrected to reflect the correct construction 4a instead of the
#' erroneous construction 3a underlying the run size in the Colbourn tables
#' (thanks to Christopher Smolen for highlighting the coverage violation).\cr
#' The cases that were checked
#' were found correct. It was not possible with the
#' available computing resources to check the actual coverage for
#' large cases, and sources apart from the Colbourn tables were neither
#' found for many of the arrays.
#'
#' The reference for Paley type constructions is Colbourn (2015),
#' and the constructions were taken from there.
#' For strength 3, only constructions
#' without doubling are included, as CAs from doubling are not competitive.
#' Small corrections to Colbourn's Table 1 were needed for primes 7,
#' 13 (strength 3 not possible) and 17 (no row can be removed).
#'
#' The reference for \code{RMScopies} constructions with \code{m=0}
#' is Raaphorst, Moura and Stevens (2014).
#' The proof that the construction is also valid for all pairs \code{(q, m)}
#' with \code{m > 0} has been provided by Mickael Bergeron using Claude.
#' It has not yet been checked independently by a human; the 3-coverage
#' for the parameter combinations in \code{rmsCOPIEScat}
#' are verified, either by exhaustive checks or in the original
#' proof by RMS 2014 (for \code{m=0}).
#' The data frame could be extended after verifying further cases.
#'
#' The creation of the objects is documented by R files in the \code{extdata}
#' folder for this package, which can be located in the directory \code{inst/extdata}
#' of the GitHub repo.
#'
#' @section Mathematical detail for the RMScopies construction:
#'
#' All arithmetic is in the Galois field GF(q), as provided by
#' \code{\link[lhs]{create_galois_field}}
#' (elements coded as 0, ..., q-1;
#' for prime q, this is arithmetic modulo q).
#'
#' \bold{Base array} (\code{m=0}, Raaphorst, Moura and Stevens 2014):
#' the k0 = q^2+q+1 columns correspond to the points \eqn{p_j}{p_j} of the
#' projective plane PG(2,q), ordered as the powers of a root of a primitive
#' cubic polynomial over GF(q) (Singer cycle).
#' The array consists of a row of zeroes, a block B1 with entries
#' \eqn{x \cdot p_j}{x.p_j} for all non-zero x in GF(q)^3, and a block B2
#' constructed in the same way with the column order reversed.
#' Three columns that are collinear in B1 are not collinear in B2, so that
#' all triples of columns are covered.
#'
#' \bold{Shifted copies} (\code{m>=1}): new column \eqn{n_i}{n_i} is a copy of
#' column \eqn{c_i = i}{c_i = i} in the zero row and in B1, and equals
#' \eqn{c_i + 1}{c_i + 1} (in GF(q)) in B2. The only uncovered 3-tuples are then
#' found in column triples \eqn{\{c_i, n_i, w\}}{{c_i, n_i, w}}, for differences
#' \eqn{n_i - c_i}{n_i - c_i} different from 0 and 1 (field elements coded 2, ..., q-1).
#'
#' \bold{Layers} (\code{m>=2}): each layer consists of q^2 runs indexed by
#' y in GF(q)^2. The q+1 directions of PG(1,q) yield q+1 mutually orthogonal
#' linear functions \eqn{y \cdot P}{y.P}. In a layer, up to q copies are active:
#' active copy i receives its own direction (\eqn{c_i = y \cdot P_s}{c_i = y.P_s},
#' \eqn{n_i = y \cdot P_s + e}{n_i = y.P_s + e}),
#' all other columns receive direction \eqn{P_0}{P_0}. Each copy is active in
#' q-2 layers, with e running through the field elements coded 2, ..., q-1. A final run (0 for the original columns,
#' 1 for the new columns) covers the image of x=0 in B2.
#' The number of layers is \eqn{L = \max(q-2, \lceil m(q-2)/q \rceil)}{L = max(q-2, ceiling(m(q-2)/q))}.
#'
#' For \code{m=1}, the q(q-2) layers' worth of runs are replaced by q^2(q-2) simple runs
#' (column \eqn{c_1=a}{c_1=a}, column \eqn{n_1 = a+e}{n_1 = a+e}, all other columns equal to s),
#' and no final run is needed.
#'
#' Run sizes: \eqn{2q^3-1}{2q^3-1} for \code{m=0},
#' \eqn{2q^3-1+q^2(q-2)}{2q^3-1+q^2(q-2)} for \code{m=1},
#' \eqn{2q^3+q^2 L}{2q^3+q^2 L} for \code{m>=2}.
#'
#' Compared to the Colbourn tables (status November 2024), this improves
#' 78 entries: v=5 (k=32,...,37), v=7 (k=58,...,65), v=8 (k=74,...,81),
#' v=9 (k=92,...,101), v=11 (k=134,...,145), v=13 (k=184,...,197),
#' v=19 (k=382,...,401); for example CA(324; 3, 32, 5) instead of 349 runs,
#' CA(350; 3, 37, 5) instead of 365 runs, and CA(2025; 3, 100, 9) instead of 2169 runs.
#'
#' @section Details:
#' For more details, see also the documentation of functions
#' \code{\link{paleyCA}}, \code{\link{cyclotomyCA}} and \code{\link{RMScopiesCA}}.
#'
#' @references Colbourn (2010), Colbourn (2015), Raaphorst, Moura and Stevens (2014)
#'
#' @examples
#' ## an impression of the content
#' head(PALEYcat)
#' tail(PALEYcat)
#' boxplot(N ~ t, PALEYcat, las=1)
#'
#' head(CYCLOTOMYcat)
#' tail(CYCLOTOMYcat)
#' xtabs(~ t + v, CYCLOTOMYcat)
#'
#' ## the following is not tested because of run times
#' ## in CRAN checks
#' \donttest{
#' fun <- function(t,k,v) eCAN(t,k,v)$CAN
#' fun <- Vectorize(fun)
#' diff <- PALEYcat[,"N"]-fun(PALEYcat[,"t"], PALEYcat[,"k"], PALEYcat[,"v"])
#' fivenum(diff)
#' ## -102.0    2.0   35.5  168.0  438.0
#'
#' ## one entry is much better than the current Colbourn table entry
#' ## coverage of this entry (t=6, k=380) has been proven
#' ## by Hriday Rungta (see https://github.com/ugroempi/CAs/issues/62).
#' ##
#' }
#'

#'@rdname DataForMathematical
"PALEYcat"

#'@rdname DataForMathematical
"CYCLOTOMYcat"

#'@rdname DataForMathematical
"rmsCOPIEScat"
