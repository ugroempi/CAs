% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ProjectionBose.R
\name{projBoseCA}
\alias{projBoseCA}
\alias{projectionBose}
\alias{N_projectionBose}
\alias{k_projectionBose}
\title{Function to create a strength 2 CA by projection of a Bose OA}
\usage{
projBoseCA(k, v, fixNA=TRUE, ...)

projectionBose(q, c=1, s=q-c, fixNA=FALSE, ...)

N_projectionBose(k, v, cmax=Inf, ...)

k_projectionBose(N, v, cmax=Inf, ...)
}
\arguments{
\item{k}{requested number of columns}

\item{v}{number of levels}

\item{fixNA}{logical; if FALSE, keeps flexible values; if TRUE, flexible values are arbitrarily fixed.}

\item{...}{currently not used}

\item{q}{a prime power}

\item{c}{an integer number, 0 <= c <= \code{q-2}; for \code{c=0}, the unchanged array created by \code{\link{SCA_Bose}} is returned.}

\item{s}{an integer number, 2 <= s <= \code{q-c}; for \code{c=0}, \code{s} is irrelevant.}

\item{N}{maximum acceptable number of runs}

\item{cmax}{the maximum number permitted for \code{c}; setting \code{cmax} prevents ridiculous values of \code{N} for cases
for which \code{k} and \code{v} are too far apart}
}
\value{
\code{projBose} returns a covering array in \code{k} columns at \code{v} levels each with \code{v} constant rows,
if this is reasonably possible by projection of a strength 2 Bose OA, which should be the case whenever
\code{v<k} and \code{k-v} not larger than seven and \code{v+(k-v-1)/2} a prime power.\cr

\code{projectionBose} returns a (\code{q^2-c} x \code{q+1+c} matrix with the first \code{q+1} columns
in \code{q-c} levels (labeled from 0 to \code{q-c-1}) and the last \code{c} columns in \code{s} levels (labeled from 0 to \code{s-1})
that is a covering array of strength 2.\cr
\code{N_projectionBose} and \code{k_projectionBose} return a named
vector with elements \code{N}, \code{k}, \code{kmax}, \code{v}, \code{q} and \code{c}.
}
\description{
A strength 2 CA with q^2-c runs with q+1 (q-c)-level columns and c s-level columns is created from a q-level Bose array arranged as SCA
}
\examples{
# create a CA(24,2,7,4)
projectionBose(5)  ## four constant rows, some flexible values
projectionBose(5, fixNA=TRUE) ## flexible values arbitrarily filled
projectionBose(5, c=2, s=2) ## 6 3-level columns (3=5-2), two 2-level columns

# create a CA(48,2,9,6)
dim(D <- projectionBose(7))
coverage(D, 2)
eCAN(2,9,6) ## two runs difference

# create a CA(63,2,9,7)
dim(projBoseCA(9,7))
eCAN(2, 9, 7)  ## four runs difference

# create a CA(166,2,15,12)
dim(projBoseCA(15,12))
eCAN(2, 15, 12)  ## optimal

# querying k or N
N_projectionBose(k=21, v=15)
k_projectionBose(N=357, v=15)
eCAN(2, 21, 15) ## the above design is the basis
eCAK(2, 357, 15) ## four more columns

}
\references{
Colbourn (2008, Theorem 2.3), Torres-Jimenez et al. (2019, Fig. 25)
}
