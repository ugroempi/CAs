% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maxconstant.R
\name{permvals}
\alias{permvals}
\alias{maxconstant}
\title{Function to make the maximum possible number of rows constant}
\usage{
maxconstant(D, verbose=0, remove=FALSE, one_is_enough=FALSE, dupcheck=FALSE, ...)
}
\arguments{
\item{D}{a uniform CA without duplicated rows, unless dupcheck is set to TRUE}

\item{...}{currently not used}

\item{verbose}{0 for neither comments nor attributes added to output objects,\cr
1 for printing messages with the row numbers used for the constant rows,\cr
2 for not printing anything but adding the list of possible row sets (that could have been made constant) as an attribute to the returned object,\cr
12 for combining both}

\item{remove}{logical, should constant rows be removed? If TRUE, the function returns a matrix with fewer rows than \code{D}.}

\item{one_is_enough}{logical; should the algorithm limit itself to making one row constant?}

\item{dupcheck}{logical, should duplicate rows be checked and removed? If FALSE, the function assumes that there are no duplicate rows.}
}
\value{
an equivalent version of \code{D} that has as many constant rows as possible,
located at the start of the array, or removed in case of \code{remove=TRUE}. \cr
For \code{verbose=2}, the attribute \code{constant_rows} is attached to the returned object;
it is a list with elements \code{design_name} and \code{row_set}, where the latter
contains the row numbers in the original design of the rows that were moved to the front.
}
\description{
applies maximum clique search for making maximally many rows constant, which
may be useful when using an array in a product or power based construction; the constant rows are moved to the front or removed (if requested)
}
\section{Details}{

The function uses a method presented by Torres-Jimenez at NIST that identifies the maximum
number of constant rows via maximum clique size in a graph, whose vertices are the rows
of \code{D} and where there are edges between rows that do not coincide in any column.
For a clique, all rows can be made constant by the equivalent transform of
swapping levels within columns. Of course, the maximum clique size is at most
the number of levels, because any further row will have to coincide with an existing one.

If the array already has the maximum possible number of constant rows,
these are moved to the first runs (or removed), and a search is avoided.\cr
If there are as many constant rows in the array as is maximally possible as asserted in a clique search,
these are again moved to the first runs (or removed).\cr
If there are constant rows that are all together in one of the largest cliques,
the first such clique is used.\cr
If there are no constant rows and only one constant row is possible, the first row is made constant.\cr
For all other cases, the last largest clique is used.\cr
The constant rows have ascending levels, starting with the lowest level.
}

\section{Warning}{

For arrays with many rows, the search for the largest cliques can take a long time.
}

\examples{
set.seed(1212) ## for array randomization
D <- lhs::createBush(4,5)
table(nvals <- lengths(lapply(1:64, function(obj) unique(D[obj,]))))
## one row is already constant:
which(nvals==1)
## 4 rows can be made constant (that is the maximum for 4 level CAs)
Dmoreconst <- maxconstant(D, verbose=1)
head(Dmoreconst)
## previous constant row 22 is now the 2nd row,
## the other clique rows were made constant and moved to positions 1, 3 and 4

## only a single row can be made constant (that is always possible)
D <- cyc(19,2) ## k=19, q=19, type 1, t=3
table(nvals <- lengths(lapply(1:19, function(obj) unique(D[obj,]))))
## no row is already constant (all have 2 distinct values)
## only a single constant row is achievable
head(maxconstant(D, verbose=1))

## two rows can be made constant
CA.22.3.12.2 <- cyc(11,2, type="3a") ## N=22, k=12, q=11, type 3a, t=3
table(nvals <- lengths(lapply(1:12, function(obj) unique(CA.22.3.12.2[obj,]))))
## no row is already constant (all have 2 distinct values)
## two constant rows are achievable
head(Dc <- maxconstant(CA.22.3.12.2, verbose=12))
## which rows of design CA.22.3.12.2 were used ?
attributes(Dc)

## example of Colbourn and Torres-Jimenez (2013)
## Figure 2 of the paper
hilf <- strsplit(c("2222001", "0022222", "1121201", "0120110",
                        "2210202", "2102122", "1200020", "0211121",
                        "1012100", "2001210", "0000001", "1111012",
                        "2222222", "2222011"), "")
plan <- do.call(rbind, lapply(hilf, as.numeric))
table(nvals <- lengths(lapply(1:14, function(obj) unique(plan[obj,]))))
## one row is already constant, two rows can be made constant
maxconstant(plan)

}
\references{
Torres-Jimenez, Acevedo-Juarez and Avila-George (2021) for the mention of the largest clique idea
}
