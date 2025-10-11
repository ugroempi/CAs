#' Creation of Covering Arrays
#' @description Creates covering arrays.
#'
#' @section Details:
#' This package constructs covering arrays, i.e., arrays that cover all \eqn{t}-ary combinations of a set of factors at least once.
#' The focus is on mathematical constructions. Initially, the package only offers arrays for which all columns have the same number of levels.
#'
#' The goal is to implement as many constructions as possible that yield arrays with small numbers of runs, as evidenced by the
#' Colbourn covering array tables, which do not provide the CAs themselves but pointers to which CA constructions yield the smallest known array
#' for which setting.
#'
#' All references of the package are listed in this file and referenced from the other documentation files.
#'
#' Within the package, available CA constructions for specific situations can be queried using
#' the function \code{\link{Ns}}, which provides sizes from available constructions, as well as the current known optimum size;
#' it is not guaranteed that all implemented constructions are incorporated, as this is not trivial for constructions with
#' intricate ingoing quantities.
#' An analogous function \code{ks} is far less advanced, as the typical application situation
#' presumably starts from a set of variables to be covered.
#'
#' Function \code{\link{bestCA}} makes use of function \code{\link{Ns}} (with the above-stated limitations)
#' and produces the best currently implemented design, including a download from internet libraries,
#' where applicable.
#'
#' Besides the construction functions, coverage properties of any array can be checked by function
#' \code{\link{coverage}} and plotted by function \code{\link{coverplot}}.
#'
#' So far, the following constructions for uniform CAs have been implemented:
#'
#' * strength 2 2-level CAs by Kleitman and Spencer (1973) and Katona (1973) in function \code{\link{KSK}},
#' * strength 3 to 6 2-level constructions via Paley matrices according to Colbourn (2015) in function \code{\link{paleyCA}},
#' * orthogonal arrays from a Bose (1938) or Bush (1952) construction in functions \code{\link{SCA_Bose}} and \code{\link{SCA_Busht}}; the latter does one more column for strength 3 with \code{q} a power of 2, according to the Sherwood, Martirosyan and Colbourn (2006) example for the usage of permutation vectors.
#' * constructions based on cyclotomy (Colbourn 2010) in function \code{\link{cyclotomyCA}} for practitioners (based on arrays listed in the Colbourn tables), and function \code{\link{cyc}} for experts who know which prime power and construction type to request for which strength.
#' * Chateauneuf and Kreher (2002) doubling in function \code{\link{CK_doublingCA}} (workhorse function \code{\link{CK_doubling}})
#' * strength 3 CAs for v=3,4,5 with up to 2v columns in 33, 88 or 185 runs (\code{\link{miscCA}})
#' * strength 3 CAs for v=q+1 with q a prime power, and k<=v, based on augmenting an ordered design with further runs (Cohen, Colbourn and Ling 2003 and 2008)
#' * the Colbourn and Torres-Jimenez (2010) power construction in function \code{\link{powerCA}}
#' * the Sherwood, Martirosyan and Colbourn (2006) construction based on SCPHFs and permutation vectors, including not only the SCPHFs from the Sherwood et al. paper but also those provided by Lanus (based on Colbourn and Lanus 2018 and Colbourn, Lanus and Sarkar 2018) in function \code{\link{scphfCA}}
#' * the construction based on CPHFs using CPHFs from Wagner, Colbourn and Simos (2022) in function \code{\link{cphfCA}}
#' * Colbourn et al. (2010, CKRS) cross-sum of an \code{N x k} code with a strength t-1 CA to yield a strength t CA
#' * recursive Bose constructions using \code{\link{productPCA}} or \code{\link{productCA}}, implemented in functions \code{\link{recursiveBose}} and \code{\link{recBoseCA}}; the variant with \code{link{productPCA}} is uniformly better than a similar construction of Hartman (2005, Theorem 7.1 and Corollary 7.2).
#' * a fusion construction based on a Bose matrix, from Colbourn (2008), in function \code{\link{fuseBoseCA}}
#' * a fusion construction based on a Bush matrix, from Colbourn (2008), in function \code{\link{fuseBushtCA}}
#' * a projection construction based on a Bose matrix, from Colbourn (2008), in functions \code{\link{projectionBose}} and \code{\link{projBoseCA}}
#' * the group-based strength 2 construction by Meagher and Stevens (2005) which uses starter vectors and a group of cycling permutations that leaves one value fixed (function \code{\link{CS_MS}})
#' * the strength 2 cover starter construction by Lobb et al. (2012) which works similarly to Meagher and Stevens (2005) but fix one or more symbols using the additive group on the non-fixed symbols (function \code{\link{CS_LCDST}})
#' * the strength 2 cover starter construction by Colbourn et al. (2006) which is simpler than Meagher and Stevens (2005) and allows a recursive product construction (function \code{\link{CS_CMMSSY}})
#' * the 2-level cover starter construction by Colbourn and Keri (2009), which creates strength 4, 5, and 6 CAs and is also closely related to the Paley type construction of Colbourn (2015)
#' * direct product construction for strength 2 CAs in function \code{\link{productCA}}, with generalizations for slightly reducing the number of runs
#' * product construction for strength 2 PCAs in \code{\link{productPCA}}, proposed by Colbourn et al. (2006), improved by Colbourn and Torres-Jimenez (2013); automated in function \code{\link{dpCA}} based on some ingredients; the construction is behind various current best CAs of the Colbourn tables, which are not yet achieved by the package.
#' * adding a column, as detailed in Theorem 3.2 of Colbourn et al. (2010, CKRS)
#' * composition of arrays, i.e., cross product of levels for arbitrary strength (e.g., Theorem 2.8 of Zhang et al. 2014), \code{\link{crossCAs}}; automated in function \code{\link{pcaCA}} based on some ingredients
#' * the CA EXtender for strength 2 3-level CAs (Torres-Jimenez, Acevedo-Juarez and Avila-George 2021), which yields the current best arrays for all scenarii (but, contrary to the 2-level case, these arase not theoretically optimal and might improve in the future) in function \code{\link{CAEX}}.
#' * identification of flexible values as in Colbourn and Torres-Jimenez (2013), and the Nayeri et al. (2013) postprocessing in function \code{\link{postopNCK}}
#'
#' So far, there are three functions for mixed level CAs: \code{MCA2} implements a strength 2
#' construction, generalized from Sherwood (2008) by Groemping (2025),
#' that is based on expanding some columns of a (mixed or)
#' uniform CA, after reducing their number of levels, by replacing the flexible values of those
#' columns with suitable small matrices and replicating the rest of those columns the corresponding
#' number of times. Its details can be found in Groemping (2025).\cr
#' \code{\link{CA_to_MCA}} takes
#' a uniform (or mixed level) array, removes levels (makes them flexible) as required.\cr
#' Both constructions may benefit from subsequent
#' removal of as many rows as possible via the Nayeri et al. (2013) construction. For
#' \code{CA_to_MCA}, this is the crucial step; for \code{MCA2}, gains are less dramatic,
#' but still often relevant.\cr
#' \code{\link{projBoseMCA}} implements the mixed level construction by Stevens, Ling and Mendelsohn (2002)
#' in the form stated by Colbourn (2008) in his Corollary 2.2, and the construction of Colbourn's Theorem 2.3.
#'
#' @section Exported objects:
#' The constructions are based on various catalogue objects and occasional arrays that can be inspected by expert users,
#' but are mainly meant to be used by the package functions:
#'
#' * \code{\link{colbournBigFrame}},
#' * \code{\link{CAEX_CAs}} and \code{\link{CAEX_lineages}},
#' * \code{\link{TJcat}} and \code{\link{TJ2level_CAs}},
#' * \code{\link{powerCTcat}},
#' * \code{\link{PCAcat}}
#' * \code{\link{DPcat}}
#' * \code{\link{SCPHFcat}}, \code{\link{SMC_SCPHFs}} and \code{\link{CL_SCPHFs}}
#' * \code{\link{CPHFcat}} and \code{\link{WCS_CPHFs}}
#' * \code{\link{MeagherStevensCombis}} and \code{\link{MeagherStevensStarters}},
#' * \code{\link{CMMSSYCombis}} and \code{\link{CMMSSYStarters}},
#' * \code{\link{LCDSTCombis}} and \code{\link{LCDSTStarters}},
#' * \code{\link{ColbournKeriCombis}} and \code{\link{ColbournKeriStarters}},
#' * \code{\link{CKRScat}} and \code{\link{CKRS_CAs}},
#' * \code{\link{WKScat}}, \code{\link{WKS_CAs}},
#' * \code{\link{PALEYcat}}, \code{\link{CYCLOTOMYcat}},
#' \code{\link{DWYERcat}}, \code{\link{NISTcat}},
#' \code{\link{miscCAcat}}, and individual arrays named in the
#' latter.
#'
#' The arrays referenced by \code{\link{DWYERcat}} and \code{\link{NISTcat}} can only be
#' accessed with an internet connection.
#'
#' The R objects are available in the \code{data} folder of the package,
#' corresponding raw data and construction code can be found in the
#' \code{extdata} folder of the package, which can be located in the
#' directory \code{inst/extdata} of the GitHub repo or, for a given
#' R installation, using \code{system.file("extdata", package="CAs")}
#' from within R.
#'
#' @author Author: Ulrike Groemping, BHT Berlin.
#'
#' @references
#' Ball, W. W. R. and Coxeter, H. S. M. (1987). Mathematical Recreations and Essays, 13th ed. New York: Dover, pp. 308-309.
#'
#' Bose, R.C. (1938). On the application of the properties of Galois fields to the problem of construction of hyper-Graeco-Latin squares. Sankhya 3, 323-338.
#'
#' Bush, K.A. (1952). Orthogonal Arrays of Index Unity. Annals of Mathematical Statistics 23, 426-434.
#'
#' Chateauneuf, M., Colbour.n, C. and Kreher, D.L. (1999). Covering Arrays of Strength Three. Des. Codes Cryptogr. 16, 235-242.
#'
#' Chateauneuf, M. and Kreher, D.L. (2002). On the state of strength‐three covering arrays. Journal of Combinatorial Designs, vol. 10, no. 4, pp. 217-238. doi: 10.1002/jcd.10002.
#'
#' Cohen, M.B., Colbourn, C.J., Ling, A.C.H. (2003). Augmenting simulated annealing to build interaction test suites, in: 14th International Symposium on Software Reliability Engineering, 2003. ISSRE 2003. IEEE, Denver, Colorado, USA, pp. 394-405. https://doi.org/10.1109/ISSRE.2003.1251061
#'
#' Cohen, M.B., Colbourn, C.J., Ling, A.C.H. (2008). Constructing strength three covering arrays with augmented annealing. Discrete Mathematics 308, 2709-2722. https://doi.org/10.1016/j.disc.2006.06.036
#'
#' Colbourn, C.J. (without year). Covering array tables: 2 ≤v ≤25, 2 ≤t≤6, t≤k ≤10000, 2005–23. \url{https://www.public.asu.edu/~ccolbou/src/tabby}.
#'
#' Colbourn, C.J., (2004). Combinatorial aspects of covering arrays. Matematiche (Catania) 59, 125-172.
#'
#' Colbourn, C.J., (2008). Strength two covering arrays: Existence tables and projection. Discrete Math., vol. 308, 772-786. https://doi.org/10.1016/j.disc.2007.07.050
#'
#' Colbourn, C. J. (2010). Covering arrays from cyclotomy, Des. Codes Cryptogr., vol. 55, no. 2, pp. 201–219. doi: 10.1007/s10623-009-9333-8.
#'
#' Colbourn, C.J. (2015). Suitable Permutations, Binary Covering Arrays, and Paley Matrices, in: Colbourn, C.J. (Ed.), Algebraic Design Theory and Hadamard Matrices, Springer Proceedings in Mathematics & Statistics. Springer International Publishing, Cham, pp. 29-42. https://doi.org/10.1007/978-3-319-17729-8_3
#'
#' Colbourn, C.J., Dougherty, R.E., Horsley, D. (2019). Distributing hash families with few rows. Theoretical Computer Science 800, 31-41. doi: 10.1016/j.tcs.2019.10.014
#'
#' Colbourn, C.J., Kéri, G., (2009). Binary Covering Arrays and Existentially Closed Graphs, in: Chee, Y.M., Li, C., Ling, S., Wang, H., Xing, C. (Eds.), Coding and Cryptology, Lecture Notes in Computer Science. Springer Berlin Heidelberg, Berlin, Heidelberg, pp. 22-33. https://doi.org/10.1007/978-3-642-01877-0_3
#'
#' Colbourn, C., Kéri, G., Rivas Soriano, P.P., Schlage-Puchta, J.-C., (2010). Covering and radius-covering arrays: constructions and classification. Discrete Appl. Math 158, 1158-1180. https://doi.org/10.1016/j.dam.2010.03.008
#'
#' Colbourn, C.J., Lanus, E., 2018. Subspace restrictions and affine composition for covering perfect hash families. Art Discrete Appl. Math. 1, #P2.03. https://doi.org/10.26493/2590-9770.1220.3a1
#'
#' Colbourn, C.J., Lanus, E., Sarkar, K., 2018. Asymptotic and constructive methods for covering perfect hash families and covering arrays. Des. Codes Cryptogr. 86, 907-937. https://doi.org/10.1007/s10623-017-0369-x
#'
#' Colbourn, C.J., Martirosyan, S.S., Mullen, G.L., Shasha, D., Sherwood, G.B., Yucas, J.L., (2006). Products of mixed covering arrays of strength two. \emph{J of Combinatorial Designs} \bold{14}, 124-138. doi: 10.1002/jcd.20065
#'
#' Colbourn, C.J., Torres-Jimenez, J. (2010). Heterogeneous hash families and covering arrays, in: Bruen, A.A., Wehlau, D.L. (Eds.), Contemporary Mathematics. American Mathematical Society, Providence, Rhode Island, pp. 3-15. https://doi.org/10.1090/conm/523/10309
#'
#' Colbourn, C.J., Torres-Jimenez, J. (2013). Profiles of covering arrays of strength two. J. Alg. Comput. 44, 31-59. doi: 10.22059/jac.2013.7914.
#'
#' Colbourn, C.J., Zhou, J. (2012). Improving Two Recursive Constructions for Covering Arrays. Journal of Statistical Theory and Practice 6, 30-47. https://doi.org/10.1080/15598608.2012.647489
#'
#' Dwyer, A. (2024). CA Database: Data base of covering arrays and related objects. \url{https://github.com/aadwyer/CA_Database}.
#'
#' Groemping, U. (2025). Generalizing a Sherwood (2008) construction for mixed level covering arrays. Report 01/2025, Reports in Mathematics, Physics and Chemistry. BHT Berlin. \url{https://www1.beuth-hochschule.de/FB_II/reports/Report-2025-001.pdf}
#'
#' Hartmann, A. (2005). Software and Hardware Testing Using Combinatorial Covering Suites. In: Golumbic, M.C., Hartman, I.BA. (eds) *Graph Theory, Combinatorics and Algorithms*. Operations Research/Computer Science Interfaces Series, vol 34. Springer-Verlag, New York.
#'
#' Ji, L. and Yin, J., 2010. Constructions of new orthogonal arrays and covering arrays of strength three. Journal of Combinatorial Theory, Series A 117, 236-247. https://doi.org/10.1016/j.jcta.2009.06.002
#'
#' Kleitman, D.J. and Spencer, J. (1973). Families of k-independent sets, Discrete Math., vol. 6, no. 3, pp. 255-262. https://doi.org/10.1016/0012-365X(73)90098-8.
#'
#' Katona, G. O. H. (1973). Two applications (for search theory and truth functions) of Sperner type theorems, Period. Math. Hung., vol. 3, no. 1, pp. 19-26. https://doi.org/10.1007/BF02018457.
#'
#' Kokkala, J. I., Meagher, K., Naserasr, R., Nurmela, K. J., Östergård, P. R. J., & Stevens, B. (2018). Classification of small strength-2 covering arrays. Zenodo. https://doi.org/10.5281/zenodo.1476059
#'
#' Lobb, J.R., Colbourn, C.J., Danziger, P., Stevens, B., Torres-Jimenez, J., (2012). Cover starters for covering arrays of strength two. Discrete Math. vol. 312, 943-956. https://doi.org/10.1016/j.disc.2011.10.026
#'
#' Martirosyan, S. and Trung, T.V. (2004). On t-Covering Arrays. Des. Codes Cryptogr., vol. 32, no. 2, pp. 323–339. https://doi.org/10.1023/B:DESI.0000029232.40302.6d.
#'
#' Meagher, K., Stevens, B., (2005). Group construction of covering arrays. J of Combinatorial Designs 13, 70-77. https://doi.org/10.1002/jcd.20035
#'
#' Meagher, K., (2005a). Group Construction of Covering Arrays: Part 2. (Unpublished technical report). Ottawa.
#'
#' Meagher, K., (2005b). Covering Arrays on Graphs: Qualitative Independence Graphs and Extremal Set Partition Theory. PhD thesis. University of Ottawa, Ottawa.
#'
#' Moura, L., Stardom, J., Stevens, B., Williams, A. (2003). Covering arrays with mixed alphabet sizes. J of Combinatorial Designs 11, 413-432. https://doi.org/10.1002/jcd.10059
#'
#' Nayeri, P., Colbourn, C.J., Konjevod, G., (2013). Randomized post-optimization of covering arrays. European Journal of Combinatorics 34, 91-103. https://doi.org/10.1016/j.ejc.2012.07.017
#'
#' NIST Covering Array Tables (last modified 2008, accessed 12 Jan 2025). \url{https://math.nist.gov/coveringarrays/}.
#'
#' Sherwood, G.B., Martirosyan, S.S., Colbourn, C.J. (2006). Covering arrays of higher strength from permutation vectors. J. of Combinatorial Designs 14, 202-213. https://doi.org/10.1002/jcd.20067
#'
#' Sherwood, G.B. (2008). Optimal and near-optimal mixed covering arrays by column expansion. Discrete Mathematics 308, 6022-6035. https://doi.org/10.1016/j.disc.2007.11.021
#'
#' Stevens, B., Ling, A. and Mendelsohn, E. (2002). A Direct Construction of Transversal Covers Using Group Divisible Designs. Ars Combinatoria 63, 145-159.
#'
#' Torres-Jimenez, J. (without year, accessed 10 Feb 2025). Covering arrays. \url{https://www.tamps.cinvestav.mx/~oc/}.
#'
#' Torres-Jimenez, J., Acevedo-Juarez, B., Avila-George, H. (2021). Covering array EXtender. Applied Mathematics and Computation 402, 126122. doi: 10.1016/j.amc.2021.126122
#'
#' Wagner, M., Colbourn, C., Simos, D.E., (2022). In-Parameter-Order strategies for covering perfect hash families. Applied Mathematics and Computation 421, 126952. https://doi.org/10.1016/j.amc.2022.126952
#'
#' Wagner, M., Kampel, L., Simos, D.E., (2021). Heuristically Enhanced IPO Algorithms for Covering Array Generation, in: Flocchini, P., Moura, L. (Eds.), Combinatorial Algorithms, Lecture Notes in Computer Science. Springer International Publishing, Cham, pp. 571–586. https://doi.org/10.1007/978-3-030-79987-8_40
#'
#' Zhang, J., Zhang, Z. and Ma, F. (2014). Automatic Generation of Combinatorial Test Data. SpringerBriefs in Computer Science. Springer Verlag, Berlin.
#'

#' @aliases 'CAs-package'
#' @importFrom grDevices rgb
#' @importFrom graphics polygon
#' @importFrom prepplot prepplot
#' @importFrom utils globalVariables
#' @importFrom sets as.set set_is_proper_subset
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
