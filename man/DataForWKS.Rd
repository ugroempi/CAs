% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataForWKS.R
\docType{data}
\name{WKScat}
\alias{WKScat}
\alias{WKS_CAs}
\title{Data for function wksCA}
\format{
\code{WKScat} is a data.frame in 43 rows with columns \code{t} (strength), \code{k} (number of CA columns),
\code{v} (number of levels) and code{N} (number of CA rows).

\code{WKS_CAs} is a list of stored strength 6 CAs for 2-level factors,
constructed by Wagner, Kampel and Simos (2021). Its element names are
the run numbers.

An object of class \code{list} of length 43.
}
\source{
The reference is Wagner, Kampel and Simos (2021), and the arrays
were downloaded at \url{https://srd.sba-research.org/data/sipo/}.\cr
The creation of the objects is documented by R files in the \code{extdata}
folder of this package, which can be located using
\code{system.file("extdata", package="CAs")}; the relatively large
raw arrays are not provided there, as they are available on the web.
}
\usage{
WKScat

WKS_CAs
}
\description{
The data.frame WKScat supports the identification of CAs from the
list WKS_CAs
}
\section{Details}{

All arrays have strength 6 and 2-level columns.
The data.frame \code{WKScat} holds information on available settings,
the list \code{WKS_CAs} holds the actual arrays.
}

\examples{
## an impression of the content
head(WKScat)
tail(WKScat)
fivenum(WKScat$N)  ## 440.0 530.5 605.0 669.0 718.0
fivenum(WKScat$k)  ## 30.0 40.5 51.0 61.5 72.0

## the following is not run because of run times
## in CRAN checks
\dontrun{
fun <- function(t,k,v) eCAN(t,k,v)$CAN
fun <- Vectorize(fun)
diff <- WKScat[,"N"]-fun(WKScat[,"t"], WKScat[,"k"], WKScat[,"v"])
fivenum(diff)
## 0 0 0 0 0
## at present, they are all optimal
}

}
\references{
Wagner, Kampel and Simos (2021)
}
\keyword{datasets}
