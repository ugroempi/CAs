#' Creation of Covering Arrays
#' @description Creates covering arrays.
#' @details This package constructs covering arrays, i.e., arrays that cover all \eqn{t}-ary combinations of a set of factors at least once.
#' The focus is on mathematical constructions. Initially, the package only offers arrays for which all columns have the same number of levels.
#'
#' The goal is to implement as many constructions as possible that yield arrays with small numbers of runs, as evidenced by the
#' Colbourn covering array tables, which do not provide the CAs themselves but pointers to which CA constructions yield the smallest known array
#' for which setting.
#'
#' All references of the package are listed in this file and referenced from the other documentation files.
#'
#' Within the package, available CA constructions for specific situations can be queried using
#' the guide functions \code{guide_CAs} (not yet implemented, might be changed).
#'
#' Besides the construction functions, coverage properties of any array can be checked by function
#' \code{\link{coverage}} and plotted by function \code{\link{coverplot}}.
#'
#' So far, the following constructions have been implemented:
#'
#' * strength 2 2-level CAs by Kleitman and Spencer (1973) and Katona (1973) in function \code{\link{KSK}},
#' * strength 3 2-level constructions via Paley-based Hadamard matrices (function \code{\link{paleyHad}})
#' * constructions based on cyclotomy (Colbourn 2010) in function \code{\link{cyc}},
#' * Chateauneuf and Kreher (2002) doubling in function \code{\link{CK_doubling}}, and also construction D from the basics part of that paper (function \code{\link{CK_constrD}})
#' * recursive Bose constructions of Hartman (2005) in function \code{\link{recursiveBose}}
#' * the group-based construction by Meagher and Stevens (2005) which uses starter vectors and a group of cycling permutations that leaves one value fixed
#' * product constructions (work in progress) in functions \code{\link{productCA1}}, \code{\link{productCA2}} and \code{\link{productPCA}}; especially the latter constructs many very good CAs (Colbourn et al. 2006).
#'
#' @author Author: Ulrike Groemping, BHT Berlin.
#'
#' @references
#' Ball, W. W. R. and Coxeter, H. S. M. (1987). Mathematical Recreations and Essays, 13th ed. New York: Dover, pp. 308-309.
#'
#' Chateauneuf, M. and Kreher, D.L. (2002). On the state of strength‐three covering arrays. Journal of Combinatorial Designs, vol. 10, no. 4, pp. 217-238. doi: 10.1002/jcd.10002.
#'
#' Colbourn, C.J. (without year). Covering array tables: 2 ≤v ≤25, 2 ≤t≤6, t≤k ≤10000, 2005–23. \url{https://www.public.asu.edu/~ccolbou/src/tabby}.
#'
#' Colbourn, C.J., Martirosyan, S.S., Mullen, G.L., Shasha, D., Sherwood, G.B., Yucas, J.L., (2006). Products of mixed covering arrays of strength two. *J of Combinatorial Designs* **14**, 124-138. doi: 10.1002/jcd.20065
#'
#' Colbourn, C. J. (2010). Covering arrays from cyclotomy, Des. Codes Cryptogr., vol. 55, no. 2, pp. 201–219. doi: 10.1007/s10623-009-9333-8.
#'
#' Hartmann, A. (2005). Software and Hardware Testing Using Combinatorial Covering Suites. In: Golumbic, M.C., Hartman, I.BA. (eds) *Graph Theory, Combinatorics and Algorithms*. Operations Research/Computer Science Interfaces Series, vol 34. Springer-Verlag, New York.
#'
#' Kleitman, D.J. and Spencer, J. (1973). Families of k-independent sets, Discrete Math., vol. 6, no. 3, pp. 255–262. doi: 10.1016/0012-365X(73)90098-8.
#'
#' Katona, G. O. H. (1973). Two applications (for search theory and truth functions) of Sperner type theorems, Period. Math. Hung., vol. 3, no. 1, pp. 19–26. doi: 10.1007/BF02018457.
#'
#' Martirosyan, S. and Trung, T.V. (2004). On t-Covering Arrays. Des. Codes Cryptogr., vol. 32, no. 2, pp. 323–339. doi: 10.1023/B:DESI.0000029232.40302.6d.
#'
#' Meagher, K., Stevens, B., (2005). Group construction of covering arrays. J of Combinatorial Designs 13, 70–77. https://doi.org/10.1002/jcd.20035
#'
#' Meagher, K., (2005a). Group Construction of Covering Arrays: Part 2. (Unpublished technical report). Ottawa.
#'
#' Meagher, K., (2005b). Covering Arrays on Graphs: Qualitative Independence Graphs and Extremal Set Partition Theory. PhD thesis. University of Ottawa, Ottawa.
#'
#' NIST Covering Array Tables (last modified 2008, accessed 12 Jan 2025). \url{https://math.nist.gov/coveringarrays/}.
#'
#' Torres-Jimenez, J. (without year, accessed 10 Feb 2025). Covering arrays. \url{https://www.tamps.cinvestav.mx/~oc/}.
#'
#' Zhang, J., Zhang, Z. and Ma, F. (2014). Automatic Generation of Combinatorial Test Data. SpringerBriefs in Computer Science. Springer Verlag, Berlin.
#'

#' @aliases 'CAs-package'
#' @importFrom grDevices rgb
#' @importFrom graphics polygon
#' @importFrom prepplot prepplot
#' @importFrom sets as.set set_is_proper_subset
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
