
<!--- DO NOT EDIT:  AUTOMATICALLY GENERATED from README.Rmd -->

# CAs

Creates covering arrays.

- **Author**: Ulrike Groemping, BHT Berlin.
- **License**: GPL-3
- **Version**: 0.2.1

## Warning

This is a very preliminary version of R package **CAs**. Changes that
break backwards compatibility can and will occur without warning.
Important changes/bug fixes are stated in the
[NEWS](https://github.com/ugroempi/CAs/blob/main/inst/NEWS.md) file, but
no guarantees!

## Installation

**CAs** is not yet on [CRAN](https://CRAN.R-project.org). If you want to
work with this very preliminary version (see above warning), you can
install the package from this repository with:

``` r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("ugroempi/CAs")
```

## Details

This package constructs covering arrays, i.e., arrays that cover all
$`t`$-ary combinations of a set of factors at least once. The focus is
on mathematical constructions. Initially, the package only offers arrays
for which all columns have the same number of levels.

The goal is to implement as many constructions as possible that yield
arrays with small numbers of runs, as evidenced by the Colbourn covering
array tables, which do not provide the CAs themselves but pointers to
which CA constructions yield the smallest known array for which setting.

All references of the package are listed in this file and referenced
from the other documentation files.

Within the package, available CA constructions for specific situations
can be queried using the guide functions <code>guide_CAs</code> (not yet
implemented, might be changed).

Besides the construction functions, coverage properties of any array can
be checked by function <code>coverage</code> and plotted by function
<code>coverplot</code>.

So far, the following constructions have been implemented:
<ul>

<li>

strength 2 2-level CAs by Kleitman and Spencer (1973) and Katona (1973)
in function <code>KSK</code>,
</li>

<li>

strength 3 to 6 2-level constructions via Paley matrices according to
Colbourn (2015) in function <code>paleyCA</code>,
</li>

<li>

constructions based on cyclotomy (Colbourn 2010) in function
<code>cyclotomyCA</code> for practitioners (based on arrays listed in
the Colbourn tables), and function <code>cyc</code> for experts who know
which prime power and construction type to request for which strength.
</li>

<li>

Chateauneuf and Kreher (2002) doubling in function
<code>CK_doubling</code>, and also construction D from the basics part
of that paper (function <code>CK_constrD</code>)
</li>

<li>

Colbourn et al. (2010, CKRS) cross-sum of an <code>N x k</code> code
with a strength t-1 CA (or is it 2 only) to yield a strength t (or is it
3 only?) CA
</li>

<li>

recursive Bose constructions of Hartman (2005) in function
<code>recursiveBose</code>
</li>

<li>

a projection construction for a Bose matrix, from Colbourn (2008), in
function <code>projectionBose</code>
</li>

<li>

the group-based strength 2 construction by Meagher and Stevens (2005)
which uses starter vectors and a group of cycling permutations that
leaves one value fixed (function <code>CS_MS</code>)
</li>

<li>

the strength 2 cover starter construction by Lobb et al. (2012) which
works similarly to Meagher and Stevens (2005) but fix one or more
symbols using the additive group on the non-fixed symbols (function
<code>CS_LCDST</code>)
</li>

<li>

direct product construction for strength 2 CAs in function
<code>productCA</code>, with generalizations for slightly reducing the
number of runs
</li>

<li>

product construction for strength 2 PCAs in <code>productPCA</code>,
proposed by Colbourn et al. (2006), improved later; this is behind
various current best CAs of the Colbourn tables.
</li>

<li>

cross product of levels for arbitrary strength (e.g., Theorem 2.8 of
Zhang et al. 2014), <code>crossCAs</code>
</li>

<li>

the CA EXtender for strength 2 3-level CAs (Torres-Jimenez,
Acevedo-Juarez and Avila-George 2021), which yields the current best
arrays for all scenarii (but, contrary to the 2-level case, these are
not theoretically optimal and might improve in the future) in function
<code>CAEX</code>.
</li>

<li>

identification of flexible values as in Colbourn and Torres-Jimenez
(2013), and the Nayeri et al. (2013) postprocessing in function
<code>postopNCK</code>
</li>

</ul>

## Exported catalogue objects

The constructions are based on various catalogue objects that can be
inspected by expert users, but are mainly meant to be used by the
package functions: <code>CAEX_CAs</code>, <code>CAEX_lineages</code>,
<code>colbournBigFrame</code>, <code>colbournCatalogue</code>,
<code>MeagherStevensCombis</code>, <code>MeagherStevensStarters</code>,
<code>LCDSTCombis</code>, <code>LCDSTStarters</code>,
<code>CKRScat</code>, <code>CKRS_CAs</code>, <code>WKScat</code>,
<code>WKS_CAs</code>, <code>DWYERcat</code>, <code>NISTcat</code>,
<code>TJcat</code>, <code>PALEYcat</code>, <code>CYCLOTOMYcat</code>.

## References

**The Colbourn Covering array tables are currently unavailable at their
usual location; temporarily, access is provided**
[here](https://github.com/ugroempi/CAs/blob/main/ColbournTables.md).

<p>

Ball, W. W. R. and Coxeter, H. S. M. (1987). Mathematical Recreations
and Essays, 13th ed. New York: Dover, pp. 308-309.
</p>

<p>

Chateauneuf, M. and Kreher, D.L. (2002). On the state of strength‐three
covering arrays. Journal of Combinatorial Designs, vol. 10, no. 4,
pp. 217-238. doi: 10.1002/jcd.10002.
</p>

<p>

Colbourn, C.J. (without year). Covering array tables: 2 ≤v ≤25, 2 ≤t≤6,
t≤k ≤10000, 2005–23.
<a href="https://www.public.asu.edu/~ccolbou/src/tabby">https://www.public.asu.edu/~ccolbou/src/tabby</a>.
</p>

<p>

Colbourn, C.J., (2008). Strength two covering arrays: Existence tables
and projection. Discrete Math., vol. 308, 772-786.
<https://doi.org/10.1016/j.disc.2007.07.050>
</p>

<p>

Colbourn, C. J. (2010). Covering arrays from cyclotomy, Des. Codes
Cryptogr., vol. 55, no. 2, pp. 201–219. doi: 10.1007/s10623-009-9333-8.
</p>

<p>

Colbourn, C.J. (2015). Suitable Permutations, Binary Covering Arrays,
and Paley Matrices, in: Colbourn, C.J. (Ed.), Algebraic Design Theory
and Hadamard Matrices, Springer Proceedings in Mathematics & Statistics.
Springer International Publishing, Cham, pp. 29–42.
<https://doi.org/10.1007/978-3-319-17729-8_3>
</p>

<p>

Colbourn, C.J., Kéri, G., (2009). Binary Covering Arrays and
Existentially Closed Graphs, in: Chee, Y.M., Li, C., Ling, S., Wang, H.,
Xing, C. (Eds.), Coding and Cryptology, Lecture Notes in Computer
Science. Springer Berlin Heidelberg, Berlin, Heidelberg, pp. 22–33.
<https://doi.org/10.1007/978-3-642-01877-0_3>
</p>

<p>

Colbourn, C., Kéri, G., Rivas Soriano, P.P., Schlage-Puchta, J.-C.,
(2010). Covering and radius-covering arrays: constructions and
classification. Discrete Appl. Math 158, 1158-1180.
<https://doi.org/10.1016/j.dam.2010.03.008>
</p>

<p>

Colbourn, C.J., Martirosyan, S.S., Mullen, G.L., Shasha, D., Sherwood,
G.B., Yucas, J.L., (2006). Products of mixed covering arrays of strength
two. <em>J of Combinatorial Designs</em> <b>14</b>, 124-138. doi:
10.1002/jcd.20065
</p>

<p>

Colbourn, C.J., Torres-Jimenez, J., (2013). Profiles of covering arrays
of strength two. J. Alg. Comput. 44, 31-59. doi: 10.22059/jac.2013.7914.
</p>

<p>

Dwyer, A., (2024). CA Database: Data base of covering arrays and related
objects.
<a href="https://github.com/aadwyer/CA_Database">https://github.com/aadwyer/CA_Database</a>.
</p>

<p>

Hartmann, A. (2005). Software and Hardware Testing Using Combinatorial
Covering Suites. In: Golumbic, M.C., Hartman, I.BA. (eds) <em>Graph
Theory, Combinatorics and Algorithms</em>. Operations Research/Computer
Science Interfaces Series, vol 34. Springer-Verlag, New York.
</p>

<p>

Kleitman, D.J. and Spencer, J. (1973). Families of k-independent sets,
Discrete Math., vol. 6, no. 3, pp. 255-262. doi:
10.1016/0012-365X(73)90098-8.
</p>

<p>

Katona, G. O. H. (1973). Two applications (for search theory and truth
functions) of Sperner type theorems, Period. Math. Hung., vol. 3, no. 1,
pp. 19–26. doi: 10.1007/BF02018457.
</p>

<p>

Lobb, J.R., Colbourn, C.J., Danziger, P., Stevens, B., Torres-Jimenez,
J., (2012). Cover starters for covering arrays of strength two. Discrete
Math. vol. 312, 943-956. <https://doi.org/10.1016/j.disc.2011.10.026>
</p>

<p>

Martirosyan, S. and Trung, T.V. (2004). On t-Covering Arrays. Des. Codes
Cryptogr., vol. 32, no. 2, pp. 323–339. doi:
10.1023/B:DESI.0000029232.40302.6d.
</p>

<p>

Meagher, K., Stevens, B., (2005). Group construction of covering arrays.
J of Combinatorial Designs 13, 70-77.
<https://doi.org/10.1002/jcd.20035>
</p>

<p>

Meagher, K., (2005a). Group Construction of Covering Arrays: Part 2.
(Unpublished technical report). Ottawa.
</p>

<p>

Meagher, K., (2005b). Covering Arrays on Graphs: Qualitative
Independence Graphs and Extremal Set Partition Theory. PhD thesis.
University of Ottawa, Ottawa.
</p>

<p>

Nayeri, P., Colbourn, C.J., Konjevod, G., (2013). Randomized
post-optimization of covering arrays. European Journal of Combinatorics
34, 91-103. <https://doi.org/10.1016/j.ejc.2012.07.017>
</p>

<p>

NIST Covering Array Tables (last modified 2008, accessed 12 Jan 2025).
<a href="https://math.nist.gov/coveringarrays/">https://math.nist.gov/coveringarrays/</a>.
</p>

<p>

Torres-Jimenez, J. (without year, accessed 10 Feb 2025). Covering
arrays.
<a href="https://www.tamps.cinvestav.mx/~oc/">https://www.tamps.cinvestav.mx/~oc/</a>.
</p>

<p>

Torres-Jimenez, J., Acevedo-Juarez, B., Avila-George, H. (2021).
Covering array EXtender. Applied Mathematics and Computation 402,
126122. doi: 10.1016/j.amc.2021.126122
</p>

<p>

Wagner, M., Kampel, L., Simos, D.E., (2021). Heuristically Enhanced IPO
Algorithms for Covering Array Generation, in: Flocchini, P., Moura, L.
(Eds.), Combinatorial Algorithms, Lecture Notes in Computer Science.
Springer International Publishing, Cham, pp. 571–586.
<https://doi.org/10.1007/978-3-030-79987-8_40>
</p>

<p>

Zhang, J., Zhang, Z. and Ma, F. (2014). Automatic Generation of
Combinatorial Test Data. SpringerBriefs in Computer Science. Springer
Verlag, Berlin.
</p>
