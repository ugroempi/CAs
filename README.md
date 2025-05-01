
<!--- DO NOT EDIT:  AUTOMATICALLY GENERATED from README.Rmd -->

# CAs

Creates covering arrays.

- **Author**: Ulrike Groemping, BHT Berlin.

- **License**: GPL-3

## Warning

This is a very preliminary version of R package **CAs**. Changes that
break backwards compatibility can and will occur without warning.

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

strength 3 2-level constructions via Paley-based Hadamard matrices
(function <code>paleyHad</code>)
</li>

<li>

constructions based on cyclotomy (Colbourn 2010) in function
<code>cyc</code>,
</li>

<li>

Chateauneuf and Kreher (2002) doubling in function
<code>CK_doubling</code>, and also construction D from the basics part
of that paper (function <code>CK_constrD</code>)
</li>

<li>

recursive Bose constructions of Hartman (2005) in function
<code>recursiveBose</code>
</li>

<li>

the group-based construction by Meagher and Stevens (2005) which uses
starter vectors and a group of cycling permutations that leaves one
value fixed
</li>

<li>

product constructions (work in progress) in functions
<code>productCA1</code>, <code>productCA2</code> and
<code>productPCA</code>; especially the latter constructs many very good
CAs (Colbourn et al. 2006).
</li>

<li>

the CA EXtender for strength 2 3-level CAs (Torres-Jimenez,
Acevedo-Juarez and Avila-George 2021), which yields the current best
arrays for all scenarii (but, contrary to the 2-level case, these are
not theoretically optimal and might improve in the future).
</li>

<li>

identification of flexible values as in Colbourn and Torres-Jimenez
(2013)
</li>

</ul>

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

Colbourn, C.J., Martirosyan, S.S., Mullen, G.L., Shasha, D., Sherwood,
G.B., Yucas, J.L., (2006). Products of mixed covering arrays of strength
two. <em>J of Combinatorial Designs</em> <strong>14</strong>, 124-138.
doi: 10.1002/jcd.20065
</p>

<p>

Colbourn, C. J. (2010). Covering arrays from cyclotomy, Des. Codes
Cryptogr., vol. 55, no. 2, pp. 201–219. doi: 10.1007/s10623-009-9333-8.
</p>

<p>

Colbourn, C.J., Torres-Jimenez, J., (2013). Profiles of covering arrays
of strength two. J. Alg. Comput. 44, 31-59. doi: 10.22059/jac.2013.7914.
</p>

<p>

Hartmann, A. (2005). Software and Hardware Testing Using Combinatorial
Covering Suites. In: Golumbic, M.C., Hartman, I.BA. (eds) <em>Graph
Theory, Combinatorics and Algorithms</em>. Operations Research/Computer
Science Interfaces Series, vol 34. Springer-Verlag, New York.
</p>

<p>

Kleitman, D.J. and Spencer, J. (1973). Families of k-independent sets,
Discrete Math., vol. 6, no. 3, pp. 255–262. doi:
10.1016/0012-365X(73)90098-8.
</p>

<p>

Katona, G. O. H. (1973). Two applications (for search theory and truth
functions) of Sperner type theorems, Period. Math. Hung., vol. 3, no. 1,
pp. 19–26. doi: 10.1007/BF02018457.
</p>

<p>

Martirosyan, S. and Trung, T.V. (2004). On t-Covering Arrays. Des. Codes
Cryptogr., vol. 32, no. 2, pp. 323–339. doi:
10.1023/B:DESI.0000029232.40302.6d.
</p>

<p>

Meagher, K., Stevens, B., (2005). Group construction of covering arrays.
J of Combinatorial Designs 13, 70–77.
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

Zhang, J., Zhang, Z. and Ma, F. (2014). Automatic Generation of
Combinatorial Test Data. SpringerBriefs in Computer Science. Springer
Verlag, Berlin.
</p>
