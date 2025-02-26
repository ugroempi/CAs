
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
$t$-ary combinations of a set of factors at least once. The focus is on
mathematical constructions. Initially, the package only offers arrays
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

So far, constructions for strength 2 2-level CAs by Kleitman and Spencer
(1973) and Katona (1973), as well as constructions based on cyclotomy
(Colbourn 2010) have been implemented.

## References

The Colbourn tables are currently unavailable at their usual location;
temporarily, access is provided
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
Colbourn, C. J. (2010). Covering arrays from cyclotomy, Des. Codes
Cryptogr., vol. 55, no. 2, pp. 201–219. doi: 10.1007/s10623-009-9333-8.
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
NIST Covering Array Tables (last modified 2008, accessed 12 Jan 2025).
<a href="https://math.nist.gov/coveringarrays/">https://math.nist.gov/coveringarrays/</a>.
</p>
<p>
Torres-Jimenez, J. (without year, accessed 10 Feb 2025). Covering
arrays.
<a href="https://www.tamps.cinvestav.mx/~oc/">https://www.tamps.cinvestav.mx/~oc/</a>.
</p>
<p>
Zhang, J., Zhang, Z. and Ma, F. (2014). Automatic Generation of
Combinatorial Test Data. SpringerBriefs in Computer Science. Springer
Verlag, Berlin.
</p>
