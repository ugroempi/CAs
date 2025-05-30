% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ProjectionBose.R
\name{projectionBose}
\alias{projectionBose}
\alias{N_projectionBose}
\alias{k_projectionBose}
\title{Function to create a strength 2 CA by projection of a Bose OA}
\usage{
projectionBose(q, fixNA=TRUE, ...)

N_projectionBose(k=NULL, v=NULL, ...)

k_projectionBose(N, v=NULL, ...)
}
\arguments{
\item{q}{a prime power}

\item{fixNA}{logical; if FALSE, keeps flexible values; if TRUE, the first \code{q}-3 rows are made constant.}

\item{...}{currently not used}

\item{k}{\code{NULL} or requested number of columns (\code{q+1}); if \code{NULL}, inferred from \code{v} as \code{k=v+2}; if specified, \code{k}-1 must be a prime power, or \code{k} will be increased to the nearest such number.}

\item{v}{\code{NULL} or number of levels; if \code{NULL}, inferred from \code{k} as \code{v=k+2} or from \code{N} as \code{sqrt(N+3)-1}; if \code{k} is also specified, it takes precedence.}

\item{N}{maximum acceptable number of runs}
}
\value{
\code{projectionBose} returns a (\code{q^2-3} x \code{q+1} matrix with columns in \code{q}-1 levels labeled from 0 to q-2 that is a covering array of strength 2.\cr
\code{N_projectionBose} and \code{k_projectionBose} return a named vector with elements \code{N}, \code{k}, \code{v} and \code{q}.
}
\description{
A strength 2 CA in q+1 (q-1)-level columns and q^2-3 runs is created
from a Bose array arranged as SCA
}
\details{
Colbourn (2008) proposed various variants of projection constructions.
This particular variant reduces the number of levels by one and the number of runs by 3,
without increasing the number of columns. It creates CA(q^2-3, 2, q+1, q-1).
The reduction in rows is achieved by using a graph coloring approach.
}
\examples{
# create a CA(22,2,6,4)
projectionBose(5)  ## two constant rows
projectionBose(5, fixNA=FALSE) ## two flexible values

# create a CA(46,2,8,6)
dim(D <- projectionBose(7))
coverage(D, 2)
eCAN(2,8,6) ## four runs difference

# create a CA(61,2,9,7)
dim(projectionBose(8))
eCAN(2, 9, 7)  ## two runs difference

# create a CA(166,2,14,12)
dim(projectionBose(13))
eCAN(2, 14, 12)  ## one run difference

# querying k or N
N_projectionBose(k=21)
eCAN(2, 24, 22) ## one run difference
eCAK(2, 526, 22) ## as good as it gets

N_projectionBose(v=6)
k_projectionBose(N=50)

k_projectionBose(N=118)
eCAN(2,12,10) ## two runs difference

}
